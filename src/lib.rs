#![feature(string_remove_matches)]

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader};
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
    chromosome: String,
    position: u32,
    id: Option<Vec<String>>,
    reference_bases: String,
    alternate_bases: Option<Vec<Option<String>>>,
    quality: f32,
    filter_status: String,
    info: Option<Vec<InfoEntry>>,
    end: Option<u32>,
    sample_info: Vec<SampleInfo>,
}

#[derive(Debug)]
pub struct SampleInfo {
    unparsed_info: String,
}

struct SampleIterator<'a> {
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
    pub fn records(&self) -> io::Result<impl Iterator<Item=io::Result<VcfRecord>> + '_> {
        let mut reader: Box<dyn BufRead> =
            if self.compressed {
                Box::new(BufReader::new(MultiGzDecoder::new(File::open(&self.path)?)))
            } else {
                Box::new(BufReader::new(File::open(&self.path)?))
            };

        let mut buf = String::new();
        loop {
            // todo we could probably do this more efficiently if we recorded the header size during parsing
            buf.clear();
            reader.read_line(&mut buf)?;

            if !buf.starts_with("##") {
                break;
            }
        }

        Ok(SampleIterator {
            reader,
            header: &self.header,
            buffer: String::with_capacity(1024),
        })
    }

    /// Parse all header and meta information in the VCF file in the reader, and return a header
    /// instance
    fn parse_header<R: BufRead>(reader: &mut R) -> io::Result<VcfHeader> {
        let mut file_version = String::with_capacity(32);

        reader.read_line(&mut file_version)?;
        assert!(file_version.starts_with("##fileformat="), "VCF file misses file format identifier");
        file_version.remove_matches("##fileformat=");

        // TODO parse rest of meta fields
        let mut buf = String::with_capacity(1024);
        loop {
            buf.clear();
            reader.read_line(&mut buf)?;

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
            file_format: file_version,
            has_end_column: end_column_present,
            sample_names: sample_column_names,
        })
    }
}

impl<'a> SampleIterator<'a> {
    pub(crate) fn parse_current_record(&self, header: &VcfHeader) -> VcfRecord {
        let mut fields = self.buffer.split('\t');

        VcfRecord {
            chromosome: fields.next().expect("VCF record empty").into(),
            position: fields.next().expect("VCF record misses POS entry").parse().expect("VCF record has malformed POS entry"),
            id: match fields.next().expect("VCF record misses ID entry") {
                "." => None,
                s => Some(s.split(';').map(|s| s.into()).collect())
            },
            reference_bases: fields.next().expect("VCF record misses REF entry").into(),
            alternate_bases: match fields.next().expect("VCF record misses ALT entry") {
                "." => None,
                s => Some(s.split(',')
                    .map(|s| match s {
                        "." => None,
                        s => Some(s.into())
                    })
                    .collect())
            },
            quality: fields.next().expect("VCF record misses QUAL entry").parse().expect("VCF record has malformed QUAL entry"),
            filter_status: fields.next().expect("VCF record misses FILTER entry").into(), // todo parse this as an enum
            info: fields.next().expect("VCF record misses INFO entry").split(';').map(|info| match info {
                "." => None,
                info => None, // todo parse
            }).collect(),
            end: if header.has_end_column { fields.next().expect(""); None } else { None }, // todo
            sample_info: if header.sample_names.is_some() {
                fields.next().expect("VCF record misses FORMAT entry");
                fields.map(|s| SampleInfo { unparsed_info: s.into() }).collect()
            } else {
                Vec::new()
            },
        }
    }
}

impl<'a> Iterator for SampleIterator<'a> {
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
            Err(e) => Some(Err(e))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn test() {
        let mut vcf_file = VcfFile::parse("run/example.vcf.gz", true).expect("failed to open VCF file");

        let start = Instant::now();
        let mut cells = 0;

        vcf_file.records().expect("failed to open VCF file")
            .take(40000)
            .for_each(|rec| {
                cells += rec.expect("failed to parse VCF record").sample_info.len()
            });
        println!("read {} cells in {:?}", cells, start.elapsed());
    }
}
