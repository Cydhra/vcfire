# Extremely fast VCF Parser

This is a work-in-progress parser for [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format) files.
The parser uses flate2 with zlib for higher performance on compressed files, 
parses many of the file contents lazily to avoid spending time on unnecessary operations,
and makes use of nightly-only API for faster string operations where applicable.

The parser does not validate most of the inputs to save runtime, and therefore panics on corrupted input files.
