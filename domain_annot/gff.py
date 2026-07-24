"""GFF3 annotation parser for gene and CDS feature coordinates."""

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, TextIO, Union
from pathlib import Path


@dataclass
class GeneFeature:
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str


def parse_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GFF3 column 9 attribute key=value string into a dictionary."""
    attributes = {}
    for item in attr_str.strip().split(";"):
        if not item or "=" not in item:
            continue
        key, val = item.split("=", 1)
        attributes[key.strip()] = val.strip()
    return attributes


def extract_gene_name(attributes: Dict[str, str]) -> Optional[str]:
    """
    Extract best gene name/identifier from GFF attributes.
    Checks Name, ID, locus_tag, protein_id in order.
    Strips common prefixes like 'gene-' or 'cds-'.
    """
    raw_name = (
        attributes.get("Name")
        or attributes.get("ID")
        or attributes.get("locus_tag")
        or attributes.get("protein_id")
    )
    if not raw_name:
        return None

    # Remove prefixes like 'gene-' or 'cds-' or 'rna-'
    clean_name = re.sub(r"^(gene-|cds-|rna-|protein-)", "", raw_name)
    return clean_name


def parse_gff(
    gff_file: Union[str, Path, TextIO],
    feature_type: str = "gene"
) -> Dict[str, GeneFeature]:
    """
    Parse a GFF3 file and return a dictionary mapping gene_name -> GeneFeature.

    Args:
        gff_file: Path to GFF3 file or open file handle.
        feature_type: Feature type to filter for (default: "gene", can also use "CDS").

    Returns:
        Dict mapping gene_name -> GeneFeature(gene_name, chrom, start, end, strand)
    """
    genes: Dict[str, GeneFeature] = {}

    if isinstance(gff_file, (str, Path)):
        handle = open(gff_file, "r")
        should_close = True
    else:
        handle = gff_file
        should_close = False

    try:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            ftype = parts[2]
            
            # Check if this line matches requested feature_type (or if requested 'gene', also accept 'CDS' if no 'gene' lines present)
            if ftype != feature_type:
                continue

            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue

            strand = parts[6]
            attr_str = parts[8]
            attrs = parse_attributes(attr_str)
            gene_name = extract_gene_name(attrs)

            if gene_name:
                genes[gene_name] = GeneFeature(
                    gene_name=gene_name,
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand
                )
    finally:
        if should_close:
            handle.close()

    return genes
