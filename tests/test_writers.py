"""Unit tests for writers module."""

from pathlib import Path
from domain_annot.domain_resolver import ResolvedDomain
from domain_annot.writers import write_tsv, write_bed, write_gff3


def test_writers(tmp_path: Path):
    d = ResolvedDomain(
        protein_id="MAB_RS20425",
        aa_start=15,
        aa_stop=87,
        interpro_acc="IPR001647",
        interpro_desc="DNA-binding HTH domain, TetR-type",
        entry_type="Domain",
        best_evalue=6.6e-11,
        go_terms="GO:0003677",
        strand="+",
        genome_start=4072070,
        genome_end=4072288,
        chrom="chr1"
    )
    domains = [d]

    # Test TSV writer
    tsv_file = tmp_path / "out.tsv"
    write_tsv(domains, tsv_file)
    tsv_content = tsv_file.read_text()
    assert "MAB_RS20425" in tsv_content
    assert "IPR001647" in tsv_content

    # Test BED writer
    bed_file = tmp_path / "out.bed"
    write_bed(domains, bed_file)
    bed_content = bed_file.read_text()
    assert "chr1\t4072069\t4072288\tMAB_RS20425:IPR001647\t0\t+" in bed_content

    # Test GFF3 writer
    gff_file = tmp_path / "out.gff3"
    write_gff3(domains, gff_file)
    gff_content = gff_file.read_text()
    assert "##gff-version 3" in gff_content
    assert "DomainAnnot\tprotein_domain\t4072070\t4072288" in gff_content
