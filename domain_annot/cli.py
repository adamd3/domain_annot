"""Command line interface for DomainAnnot."""

import sys
from pathlib import Path
from typing import Optional, List

import click

from domain_annot import __version__
from domain_annot.gff import parse_gff
from domain_annot.interpro import parse_interpro_entry_list, parse_interpro_tsv
from domain_annot.ncbi import fetch_genbank_by_accession, fetch_interpro_entry_list, parse_genbank, export_fasta
from domain_annot.domain_resolver import resolve_domains
from domain_annot.writers import write_tsv, write_bed, write_gff3


@click.group()
@click.version_option(version=__version__, prog_name="domain-annot")
def main():
    """DomainAnnot: Bacterial Protein Domain Annotation based on Genome Sequence."""
    pass


@main.command(name="fetch")
@click.option("-a", "--accession", required=True, help="NCBI accession ID (e.g. NC_010397.1)")
@click.option("-o", "--outdir", default="data", show_default=True, help="Output directory for downloaded files")
def fetch_cmd(accession: str, outdir: str):
    """Download GenBank annotation file from NCBI and EBI InterPro entry list."""
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    click.echo(f"[*] Downloading GenBank record for accession '{accession}'...")
    gb_file = out_path / f"{accession}.gbk"
    fetch_genbank_by_accession(accession, gb_file)
    click.echo(f"[+] GenBank record saved to: {gb_file}")

    click.echo("[*] Extracting protein FASTA from GenBank record...")
    res = parse_genbank(gb_file)
    fasta_file = out_path / f"{accession}_proteins.nostop.fasta"
    export_fasta(res.protein_fasta, fasta_file)
    click.echo(f"[+] Extracted {len(res.protein_fasta)} proteins to: {fasta_file}")

    click.echo("[*] Downloading EBI InterPro entry list...")
    entry_file = out_path / "interpro_entry_list.txt"
    fetch_interpro_entry_list(entry_file)
    click.echo(f"[+] InterPro entry list saved to: {entry_file}")

    click.echo("[✓] Data fetching complete.")


@main.command(name="process")
@click.option("-i", "--interpro-tsv", required=True, type=click.Path(exists=True), help="Path to InterProScan TSV output")
@click.option("-e", "--entry-list", type=click.Path(), help="Path to interpro_entry_list.txt (will fetch if not provided)")
@click.option("-g", "--genbank", type=click.Path(exists=True), help="Path to GenBank file for genomic mapping")
@click.option("--gff", type=click.Path(exists=True), help="Path to GFF3 file for genomic mapping (if GenBank not provided)")
@click.option("-o", "--outdir", default="results", show_default=True, help="Output directory for domain annotations")
@click.option("-prefix", "--prefix", default="domains", show_default=True, help="Output filename prefix")
def process_cmd(
    interpro_tsv: str,
    entry_list: Optional[str],
    genbank: Optional[str],
    gff: Optional[str],
    outdir: str,
    prefix: str
):
    """Process InterProScan results and generate genomic domain mapping files."""
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    # 1. Handle InterPro Entry List
    if entry_list and Path(entry_list).exists():
        entry_file = Path(entry_list)
    else:
        entry_file = out_path / "interpro_entry_list.txt"
        if not entry_file.exists():
            click.echo("[*] Downloading EBI InterPro entry list...")
            fetch_interpro_entry_list(entry_file)
    
    entries = parse_interpro_entry_list(entry_file)
    click.echo(f"[*] Loaded {len(entries)} InterPro entry classifications.")

    # 2. Parse InterProScan TSV
    click.echo(f"[*] Parsing InterProScan results from: {interpro_tsv}")
    hits = parse_interpro_tsv(interpro_tsv)
    click.echo(f"[*] Loaded {len(hits)} raw InterProScan hits.")

    # 3. Load Genomic Coordinate Mapping (from GenBank or GFF)
    genes = None
    if genbank:
        click.echo(f"[*] Extracting genomic gene coordinates from GenBank: {genbank}")
        gb_res = parse_genbank(genbank)
        genes = gb_res.genes
        click.echo(f"[*] Loaded coordinates for {len(genes)} genes.")
    elif gff:
        click.echo(f"[*] Parsing genomic gene coordinates from GFF3: {gff}")
        genes = parse_gff(gff)
        click.echo(f"[*] Loaded coordinates for {len(genes)} genes.")

    # 4. Resolve Domain Boundaries & Specificity Upgrades
    click.echo("[*] Resolving domain boundaries and specificity upgrades...")
    resolved = resolve_domains(hits, entries, genes=genes)
    click.echo(f"[+] Resolved {len(resolved)} non-redundant domain annotations across proteins.")

    # 5. Export Writers
    tsv_out = out_path / f"{prefix}.tsv"
    write_tsv(resolved, tsv_out)
    click.echo(f"[+] Saved TSV domain annotations: {tsv_out}")

    if genes:
        bed_out = out_path / f"{prefix}.bed"
        write_bed(resolved, bed_out)
        click.echo(f"[+] Saved BED domain annotations: {bed_out}")

        gff_out = out_path / f"{prefix}.gff3"
        write_gff3(resolved, gff_out)
        click.echo(f"[+] Saved GFF3 domain annotations: {gff_out}")

    click.echo("[✓] Processing complete!")


if __name__ == "__main__":
    main()
