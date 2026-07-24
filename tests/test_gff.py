"""Unit tests for GFF3 parser."""

import io
from domain_annot.gff import parse_gff, parse_attributes, extract_gene_name


SAMPLE_GFF = """##gff-version 3
chr1\tRefSeq\tgene\t1\t1476\t.\t+\t.\tID=gene-MAB_RS00150;Name=MAB_RS00150;locus_tag=MAB_RS00150
chr1\tRefSeq\tgene\t1553\t1747\t.\t-\t.\tID=gene-MAB_RS00155;Name=MAB_RS00155
chr1\tRefSeq\tCDS\t1\t1476\t.\t+\t0\tID=cds-MAB_RS00150;parent=gene-MAB_RS00150
"""


def test_parse_attributes():
    attrs = parse_attributes("ID=gene-MAB_RS00150;Name=MAB_RS00150;locus_tag=MAB_RS00150")
    assert attrs["ID"] == "gene-MAB_RS00150"
    assert attrs["Name"] == "MAB_RS00150"
    assert attrs["locus_tag"] == "MAB_RS00150"


def test_extract_gene_name():
    attrs = {"ID": "gene-MAB_RS00150", "Name": "MAB_RS00150"}
    assert extract_gene_name(attrs) == "MAB_RS00150"

    attrs_cds = {"ID": "cds-MAB_RS00150"}
    assert extract_gene_name(attrs_cds) == "MAB_RS00150"


def test_parse_gff():
    gff_handle = io.StringIO(SAMPLE_GFF)
    genes = parse_gff(gff_handle, feature_type="gene")

    assert len(genes) == 2
    assert "MAB_RS00150" in genes
    assert genes["MAB_RS00150"].start == 1
    assert genes["MAB_RS00150"].end == 1476
    assert genes["MAB_RS00150"].strand == "+"

    assert "MAB_RS00155" in genes
    assert genes["MAB_RS00155"].strand == "-"
    assert genes["MAB_RS00155"].start == 1553
    assert genes["MAB_RS00155"].end == 1747
