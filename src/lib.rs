#![feature(string_remove_matches)]
#![feature(slice_internals)]

extern crate core;

use core::slice::memchr::memchr;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read};

use flate2::read::MultiGzDecoder;

pub struct VcfFile {
    path: String,
    compressed: bool,
    pub header: VcfHeader,
}

pub struct VcfHeader {
    pub file_format: String,
    pub has_end_column: bool,
    pub sample_names: Option<Vec<String>>,

    // size of the entire header in bytes
    size: usize,
}

#[derive(Debug)]
pub enum InfoEntry {
    AncestralAllele(String),
    AlleleCount(Vec<u32>),
    TotalAlleleReadDepth(Vec<u32>),
    ForwardAlleleReadDepth(Vec<u32>),
    ReverseAlleleReadDepth(Vec<u32>),
    AlleleFrequency(Vec<f32>),
    AlleleNumber(u32),
    RmsBaseQuality(f32),
    Cigar(Vec<String>),
    SNPDatabaseMembership,
    CombinedDepth(u32),
    End(u32),
    HapMap2,
    HapMap3,
    RmsMappingQuality(f32),
    MapQReads(u32),
    SamplesWithData(u32),
    StrandBias(u32, u32, u32, u32),
    Somatic,
    Validated,
    Flag1000G,
    NonStandard(NonStandardInfoValue),
}

#[derive(Debug)]
pub enum NonStandardInfoValue {
    NoValue,
    SingleValue(String),
    ValueList(Vec<String>),
}

#[derive(Debug)]
pub struct VcfRecord {
    pub chromosome: Option<String>,
    pub position: Option<u32>,
    pub id: Option<Vec<String>>,
    pub reference_bases: Option<String>,
    pub alternate_bases: Option<Vec<Option<String>>>,
    pub quality: Option<f32>,
    pub filter_status: Option<String>,
    pub info: Option<Vec<Option<InfoEntry>>>,
    pub end: Option<u32>,
    pub sample_info: Option<SampleInfo>,
}

#[derive(Debug)]
pub struct SampleInfo {
    pub format: Vec<String>,
    unparsed_info: String,
}

#[derive(Debug)]
pub struct Sample<'a> {
    unparsed_info: &'a str,
}

pub mod parser_capabilities {
    pub const CHROM: u16 = 0b0000000000000001;
    pub const POS: u16 = 0b0000000000000010;
    pub const ID: u16 = 0b0000000000000100;
    pub const REF: u16 = 0b0000000000001000;
    pub const ALT: u16 = 0b0000000000010000;
    pub const QUAL: u16 = 0b0000000000100000;
    pub const FILTER: u16 = 0b0000000001000000;
    pub const INFO: u16 = 0b0000000010000000;
    pub const END: u16 = 0b0000000100000000;
    pub const SAMPLES: u16 = 0b0000001000000000;

    pub const ALL: u16 = 0b1111111111;
}

struct SampleIterator<'a, const FLAGS: u16> {
    reader: Box<dyn BufRead>,
    header: &'a VcfHeader,
    buffer: String,
}

impl VcfFile {
    /// Parse the header of a VCF file. The file handle will be closed after the header is parsed.
    /// Accessing records will open new file handles.
    pub fn parse(path: &str, compressed: bool) -> io::Result<VcfFile> {
        let header = if compressed {
            Self::parse_header(&mut BufReader::new(MultiGzDecoder::new(File::open(path)?)))?
        } else {
            Self::parse_header(&mut BufReader::new(File::open(path)?))?
        };

        Ok(VcfFile {
            path: String::from(path),
            compressed,
            header,
        })
    }

    // Open the VCF file and get a sequential lazy iterator over all samples
    pub fn records<const FLAGS: u16>(
        &self,
    ) -> io::Result<impl Iterator<Item = io::Result<VcfRecord>> + '_> {
        let mut reader: Box<dyn BufRead> = if self.compressed {
            Box::new(BufReader::new(MultiGzDecoder::new(File::open(&self.path)?)))
        } else {
            Box::new(BufReader::new(File::open(&self.path)?))
        };

        // we need to skip the header, and sadly GzReader doesnt provide functionality to read without
        // allocation
        let mut buf = vec![0; self.header.size];
        reader.read_exact(&mut buf)?;

        Ok(SampleIterator::<FLAGS> {
            reader,
            header: &self.header,
            buffer: String::with_capacity(1024),
        })
    }

    /// Parse all header and meta information in the VCF file in the reader, and return a header
    /// instance
    fn parse_header<R: BufRead>(reader: &mut R) -> io::Result<VcfHeader> {
        let mut file_version = String::with_capacity(32);
        let mut header_size = 0;

        header_size += reader.read_line(&mut file_version)?;
        assert!(
            file_version.starts_with("##fileformat="),
            "VCF file misses file format identifier"
        );
        file_version.remove_matches("##fileformat=");

        // trim newline
        debug_assert!(file_version.pop() == Some('\n'));
        if file_version.ends_with('\r') {
            file_version.pop();
        }

        // TODO parse rest of meta fields
        let mut buf = String::with_capacity(1024);
        loop {
            buf.clear();
            header_size += reader.read_line(&mut buf)?;

            if !buf.starts_with("##") {
                break;
            }
        }

        // parse header line
        assert!(buf.starts_with("#"), "VCF file misses header line");
        let mut head = buf.split('\t');

        let mut optional_column = head.nth(8);
        let mut end_column_present = false;
        let mut sample_column_names = None;

        if let Some(column_name) = optional_column {
            if column_name.starts_with("E") {
                end_column_present = true;

                // check next column
                optional_column = head.next();
            }
        }

        if let Some(column_name) = optional_column {
            if column_name.starts_with("F") {
                sample_column_names = Some(head.map(|s| String::from(s)).collect::<Vec<String>>())
            }
        }

        Ok(VcfHeader {
            size: header_size,
            file_format: file_version,
            has_end_column: end_column_present,
            sample_names: sample_column_names,
        })
    }
}

impl<'a, const FLAGS: u16> SampleIterator<'a, FLAGS> {
    pub(crate) fn parse_current_record(&self, header: &VcfHeader) -> VcfRecord {
        let fields_without_samples =
            8 + header.has_end_column as usize + header.sample_names.is_some() as usize;

        let mut fields = self.buffer.splitn(fields_without_samples + 1, '\t');

        VcfRecord {
            chromosome: if FLAGS & parser_capabilities::CHROM != 0 {
                Some(fields.next().expect("VCF record empty").into())
            } else {
                fields.next();
                None
            },
            position: if FLAGS & parser_capabilities::POS != 0 {
                Some(
                    fields
                        .next()
                        .expect("VCF record misses POS entry")
                        .parse()
                        .expect("VCF record has malformed POS entry"),
                )
            } else {
                fields.next();
                None
            },
            id: if FLAGS & parser_capabilities::ID != 0 {
                match fields.next().expect("VCF record misses ID entry") {
                    "." => None,
                    s => Some(s.split(';').map(|s| s.into()).collect()),
                }
            } else {
                fields.next();
                None
            },
            reference_bases: if FLAGS & parser_capabilities::REF != 0 {
                Some(fields.next().expect("VCF record misses REF entry").into())
            } else {
                fields.next();
                None
            },
            alternate_bases: if FLAGS & parser_capabilities::ALT != 0 {
                match fields.next().expect("VCF record misses ALT entry") {
                    "." => None,
                    s => Some(
                        s.split(',')
                            .map(|s| match s {
                                "." => None,
                                s => Some(s.into()),
                            })
                            .collect(),
                    ),
                }
            } else {
                fields.next();
                None
            },
            quality: if FLAGS & parser_capabilities::QUAL != 0 {
                Some(
                    fields
                        .next()
                        .expect("VCF record misses QUAL entry")
                        .parse()
                        .expect("VCF record has malformed QUAL entry"),
                )
            } else {
                fields.next();
                None
            },
            filter_status: if FLAGS & parser_capabilities::FILTER != 0 {
                Some(
                    fields
                        .next()
                        .expect("VCF record misses FILTER entry")
                        .into(),
                ) // todo parse this as an enum
            } else {
                fields.next();
                None
            },
            info: if FLAGS & parser_capabilities::INFO != 0 {
                Some(
                    fields
                        .next()
                        .expect("VCF record misses INFO entry")
                        .split(';')
                        .map(|info| match info {
                            "." => None,
                            info => None, // todo parse
                        })
                        .collect(),
                )
            } else {
                fields.next();
                None
            },
            end: if FLAGS & parser_capabilities::END != 0 {
                if header.has_end_column {
                    fields.next().expect("");
                    None // todo parse
                } else {
                    None
                }
            } else {
                if header.has_end_column {
                    fields.next();
                }
                None
            },
            sample_info: if FLAGS & parser_capabilities::SAMPLES != 0 {
                if header.sample_names.is_some() {
                    Some(SampleInfo {
                        format: fields.next().expect("VCF record misses FORMAT entry").split(':').map(|s| s.into()).collect(),
                        unparsed_info: fields
                            .next()
                            .expect("VCF record misses sample info entries")
                            .trim()
                            .into(),
                    })
                } else {
                    None
                }
            } else {
                None
            },
        }
    }
}

impl<'a, const FLAGS: u16> Iterator for SampleIterator<'a, FLAGS> {
    type Item = io::Result<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.buffer.clear();
        let result = self.reader.read_line(&mut self.buffer);
        match result {
            Ok(size) => {
                if size == 0 {
                    return None;
                }

                Some(Ok(self.parse_current_record(&self.header)))
            }
            Err(e) => Some(Err(e)),
        }
    }
}

impl SampleInfo {
    pub fn samples(&self) -> impl Iterator<Item=Sample<'_>> {
        fast_split(&self.unparsed_info, '\t' as u8)
            .map(|s| Self::parse_sample(s))
    }

    fn parse_sample(text: &str) -> Sample<'_> {
        Sample {
            unparsed_info: text,
        }
    }
}

impl<'a> Sample<'a> {
    pub fn entries(&self) -> impl Iterator<Item=&'_ str> {
        self.unparsed_info.split(':')
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;

    #[test]
    fn test() {
        let vcf_file = VcfFile::parse("run/example.vcf.gz", true).expect("failed to open VCF file");

        let start = Instant::now();
        let mut cells = 0;

        vcf_file
            .records::<{ parser_capabilities::SAMPLES }>()
            .expect("failed to open VCF file")
            .take(40000)
            .for_each(|rec| {
                cells += rec
                    .expect("failed to parse VCF record")
                    .sample_info
                    .expect("")
                    .unparsed_info
                    .split('\t')
                    .count()
            });
        println!("read {} cells in {:?}", cells, start.elapsed());
    }
}

struct FastSplitIter<'a> {
    text: &'a str,
    delim: u8,
    start: usize,
}

impl<'a> Iterator for FastSplitIter<'a> {
    type Item = &'a str;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start >= self.text.len() {
            return None;
        }

        let end = memchr(self.delim, &self.text.as_bytes()[self.start..]).unwrap_or(self.text.len() - self.start);
        let result = &self.text[self.start..self.start + end];
        self.start += end + 1;
        Some(result)
    }
}

pub fn fast_split(text: &str, delim: u8) -> impl Iterator<Item=&str> {
    FastSplitIter {
        text,
        delim,
        start: 0,
    }
}
