"""NCBI GenBank fetching and parsing module."""

import os
import re
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, TextIO
from dataclasses import dataclass

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from domain_annot.gff import GeneFeature

EBI_ENTRY_LIST_URL = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list"


@dataclass
class GenBankExtractionResult:
    genes: Dict[str, GeneFeature]           # protein_id / locus_tag -> GeneFeature
    protein_fasta: Dict[str, str]           # protein_id / locus_tag -> amino acid sequence (nostop)


def parse_genbank(
    gb_file: Union[str, Path, TextIO]
) -> GenBankExtractionResult:
    """
    Parse a GenBank file (.gb/.gbk) to extract gene coordinates and protein sequences.

    Args:
        gb_file: Path to GenBank file or open file handle.

    Returns:
        GenBankExtractionResult containing gene coordinates mapping and protein FASTA dict.
    """
    genes: Dict[str, GeneFeature] = {}
    protein_fasta: Dict[str, str] = {}

    if isinstance(gb_file, (str, Path)):
        handle = open(gb_file, "r")
        should_close = True
    else:
        handle = gb_file
        should_close = False

    try:
        records = SeqIO.parse(handle, "genbank")
        for record in records:
            chrom = record.id
            for feature in record.features:
                if feature.type == "CDS":
                    qualifiers = feature.qualifiers
                    
                    # Extract identifier: locus_tag or protein_id or gene
                    protein_id = None
                    if "locus_tag" in qualifiers:
                        protein_id = qualifiers["locus_tag"][0]
                    elif "protein_id" in qualifiers:
                        protein_id = qualifiers["protein_id"][0]
                    elif "gene" in qualifiers:
                        protein_id = qualifiers["gene"][0]

                    if not protein_id:
                        continue

                    # Extract coordinates (convert from BioPython 0-based to 1-based inclusive)
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "+"

                    genes[protein_id] = GeneFeature(
                        gene_name=protein_id,
                        chrom=chrom,
                        start=start,
                        end=end,
                        strand=strand
                    )

                    # Extract translation if available
                    if "translation" in qualifiers:
                        seq = qualifiers["translation"][0].strip()
                        # Ensure trailing stop asterisk is removed if present
                        if seq.endswith("*"):
                            seq = seq[:-1]
                        protein_fasta[protein_id] = seq

    finally:
        if should_close:
            handle.close()

    return GenBankExtractionResult(genes=genes, protein_fasta=protein_fasta)


def export_fasta(protein_fasta: Dict[str, str], out_path: Union[str, Path]) -> Path:
    """
    Write protein FASTA dictionary to a file.

    Args:
        protein_fasta: Dict mapping header_id -> amino_acid_sequence.
        out_path: File path to save FASTA.

    Returns:
        Path to written FASTA file.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for seq_id, seq in protein_fasta.items():
            f.write(f">{seq_id}\n{seq}\n")
    return out_path


def fetch_interpro_entry_list(out_path: Union[str, Path]) -> Path:
    """
    Download interpro_entry_list.txt from EBI FTP server.

    Args:
        out_path: Output file path.

    Returns:
        Path to downloaded file.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(EBI_ENTRY_LIST_URL, out_path)
    return out_path


def fetch_genbank_by_accession(accession: str, out_path: Union[str, Path]) -> Path:
    """
    Download a GenBank record from NCBI Entrez by accession ID (e.g. NC_010397.1 or CP000656.1).

    Args:
        accession: NCBI Accession ID.
        out_path: Output file path.

    Returns:
        Path to downloaded GenBank file.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=gbwithparts&retmode=text"
    urllib.request.urlretrieve(url, out_path)
    return out_path
