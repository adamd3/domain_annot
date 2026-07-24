"""Unit tests for NCBI GenBank parser and helper functions."""

import io
from pathlib import Path
from domain_annot.ncbi import parse_genbank, export_fasta


SAMPLE_GENBANK = """LOCUS       NC_000913               1000 bp    DNA     circular BCT 24-JUL-2026
DEFINITION  Escherichia coli str. K-12 substr. MG1655, complete genome.
ACCESSION   NC_000913
VERSION     NC_000913.3
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Escherichia coli str. K-12 substr. MG1655"
     CDS             190..255
                     /locus_tag="b0001"
                     /gene="thrL"
                     /protein_id="NP_414542.1"
                     /translation="MKRISTTITTTITITTGNGAG*"
     CDS             complement(337..800)
                     /locus_tag="b0002"
                     /gene="thrA"
                     /protein_id="NP_414543.1"
                     /translation="MRVLKFGGTSVANAERFLRVADILE*"
ORIGIN
        1 agcttttatt ctgactgcaa cgggcaatat gtctctgtgt ggattaaaaa aagagtgtct
//
"""


def test_parse_genbank():
    gb_handle = io.StringIO(SAMPLE_GENBANK)
    res = parse_genbank(gb_handle)

    assert len(res.genes) == 2
    assert len(res.protein_fasta) == 2

    # Check positive strand gene b0001
    assert "b0001" in res.genes
    assert res.genes["b0001"].chrom == "NC_000913.3"
    assert res.genes["b0001"].start == 190
    assert res.genes["b0001"].end == 255
    assert res.genes["b0001"].strand == "+"
    # Ensure trailing stop asterisk is removed
    assert res.protein_fasta["b0001"] == "MKRISTTITTTITITTGNGAG"

    # Check negative strand gene b0002
    assert "b0002" in res.genes
    assert res.genes["b0002"].strand == "-"
    assert res.genes["b0002"].start == 337
    assert res.genes["b0002"].end == 800
    assert res.protein_fasta["b0002"] == "MRVLKFGGTSVANAERFLRVADILE"


def test_export_fasta(tmp_path: Path):
    fasta_data = {"geneA": "MVK", "geneB": "MARS"}
    out_file = tmp_path / "test.fasta"
    export_fasta(fasta_data, out_file)

    content = out_file.read_text()
    assert ">geneA\nMVK\n" in content
    assert ">geneB\nMARS\n" in content
