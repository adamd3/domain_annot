# DomainAnnot: Bacterial Protein Domain Annotation based on Genome Sequence

`domain_annot` is a Python package and CLI application that builds a genomic map of functional domains in bacterial proteins. It takes protein domain predictions (from InterProScan), resolves non-redundant domain boundaries with specificity upgrades, and translates domain amino acid positions into genomic nucleotide coordinates.

---

## Features

- **GenBank & NCBI Integration**: Download reference genome annotations directly via NCBI Accession/TaxID and extract both protein FASTA (with stop codons stripped) and genomic coordinates.
- **InterPro Boundary Resolution**: Group raw InterProScan hits, collapse redundant boundaries, and upgrade domain specificity based on official EBI entry type priorities (`Family` > `Domain` > `Homologous_superfamily`).
- **Overlap Reduction**: Automatically filter out redundant domain predictions exceeding a 50% overlap threshold.
- **Strand-Aware Coordinate Mapping**: Accurately translate amino acid positions into 1-based genomic nucleotide coordinates for both `+` and `-` strand genes.
- **Multi-Format Exports**: Output results as clean TSV, 6-column BED, or standard GFF3 files for genome browsers.

---

## Installation

### From Source
```bash
git clone https://github.com/your-username/domain_annot.git
cd domain_annot

# Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install domain_annot
pip install -e .
```

---

## Quickstart & Usage

### 1. Download GenBank & InterPro Data (`fetch`)
Download a complete GenBank record and extract the protein sequence FASTA along with the EBI InterPro entry list:
```bash
domain-annot fetch -a NC_010397.1 -o data/
```

### 2. Process InterProScan Results (`process`)
Process InterProScan results (`.tsv`) and map resolved domains onto genomic coordinates using a GenBank record:
```bash
domain-annot process \
  -i path/to/interproscan_results.tsv \
  -g data/NC_010397.1.gbk \
  -o results/ \
  --prefix my_genome_domains
```

Or process using a separate GFF3 annotation file:
```bash
domain-annot process \
  -i path/to/interproscan_results.tsv \
  --gff path/to/annotations.gff3 \
  -o results/
```

---

## Python API Usage

You can also use `domain_annot` directly in your Python code:

```python
from domain_annot.ncbi import parse_genbank, parse_interpro_entry_list, parse_interpro_tsv
from domain_annot.domain_resolver import resolve_domains
from domain_annot.writers import write_tsv

# Load GenBank record for coordinate mapping
gb_result = parse_genbank("data/genome.gbk")

# Parse InterPro entries and InterProScan hits
entry_list = parse_interpro_entry_list("data/interpro_entry_list.txt")
hits = parse_interpro_tsv("data/interpro_results.tsv")

# Resolve domain boundaries and map coordinates
resolved = resolve_domains(hits, entry_list, genes=gb_result.genes)

# Export to TSV
write_tsv(resolved, "results/domains.tsv")
```

---

## Testing

Run the `pytest` test suite:
```bash
pytest -v
```

---

## License

MIT License
