# PyMS

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub tag](https://img.shields.io/github/tag/aretaon/PyMS.svg)](https://GitHub.com/aretaon/PyMS/tags/)
[![Github all releases](https://img.shields.io/github/downloads/aretaon/PyMS/total.svg)](ttps://GitHub.com/Naereen/StrapDown.js/releases)
[![GitHub license](https://img.shields.io/github/license/aretaon/PyMS.svg)](https://github.com/aretaon/PyMs/master/LICENSE)


This is a bunch of python scripts generated during biochemistry and/or structural proteomics and/or mass spectrometry daily work in general.
All files are work in progress and no liability is taken for them to work as intended.
Files are generally either Python modules that can be used for scripting or more general purpose programmes that come with a command line argument parser and usually also with a compiled exe for use on Windows PC with Python not available.

### Shell/Command Line Tools (Located at ./scripts)

#### waterMetadataPlotter

Raw files generated from Waters instruments come as a folder labeled FILE.raw containing all the spectra information including metadata on which parameters the measurement was started with, descriptions and so on.
In case of multiple raw files, keeping track of which file contains which measurement can be tedious. This tool summarises all raw files present in a folder by generating a single csv file from their metadata.

```
usage: watersMetadataPlotter.exe [-h] [-d DIR] [-o OUTPATH] [-a]

optional arguments:
  -h, --help            show this help message and exit
  -d DIR, --dir DIR     Path to the directory containing the raw files
  -o OUTPATH, --outpath OUTPATH
                        Path to a folder to write the output
  -a, --all             Print all columns and not only those that differ
                        between the samples
```
The script generates a csv file named FOLDERNAME_complete_tune.csv that can be viewed and edited in any spreadsheet software.

#### plotSequenceCoverageIntensity

Generates a scatter plot of the sum of MaxQuant peptide intensities (from evidence.txt) vs the position of each amino acid in the sequence. Thereby sequence coverage and for example presence degradation products can easily be investigated.

```
usage: plotSequenceCoverageIntensity.exe [-h] txt_path target_name fasta_path

positional arguments:
  txt_path     Path to the MaxQuant txt folder
  target_name  Name of the target protein as stored in the fasta file
  fasta_path   Path to fasta file

optional arguments:
  -h, --help   show this help message and exit
```
### Modules intended to be loaded into python scripts (Located at ./modules)
