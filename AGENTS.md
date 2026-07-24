# DomainAnnot: Bacterial Protein Domain Annotation based on Genome Sequence

DomainAnnot is an application that builds a genomic map of the locations of
functional domains in bacterial proteins. It uses
[InterProScan](https://www.ebi.ac.uk/interpro/) to scan protein sequences for
regions of functional homology, then assigns functional categories to
the corresponding genomic regions in which they are encoded.

## Architecture

The shell script (run_interproscan.sh) shows how the bash-based pipeline was
run. The R script was used to parse the results for an example run.
Our objective is to build a Python package and CLI tool (`domain_annot`)
(replacing the steps carried out in the R script with modular Python code).
We will convert the steps performed in the R script into clean, modular Python modules
and provide a command-line interface and Python library API; possibly
adding additional features and expanded capabilities.

## Development Standards

- Write clean, modular functions with type hints.
- Maintain snake_case naming conventions for functions and variables.
- Keep dependencies minimal; prefer standard libraries where practical.

# Testing & Verification Loop

- Unit tests live in the `tests/` directory.
- Always run `pytest` via terminal before declaring a task complete.
- Do not modify existing API/function contracts without asking.

## TODO / plans

- Build out `domain_annot` Python package & CLI.
- Set up unit testing using `pytest`.
- Download the Genbank-format annotation from NCBI based on taxonomy ID / Accession
  supplied by the user via CLI parameter or API.
  Extract both the protein sequences and GFF3 annotation from the Genbank file,
  and fetch `interpro_entry_list.txt`.
- Create git repo and GitHub repository (add README etc).

