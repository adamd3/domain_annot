"""Translate amino acid domain coordinates to genomic nucleotide coordinates."""

from typing import Tuple, Optional


def map_aa_to_genome_coords(
    gene_start: int,
    gene_end: int,
    strand: str,
    aa_start: int,
    aa_stop: int
) -> Tuple[Optional[int], Optional[int]]:
    """
    Translate 1-based inclusive amino acid domain coordinates (aa_start, aa_stop)
    into 1-based genomic nucleotide coordinates (genome_start, genome_end).

    Args:
        gene_start: 1-based start coordinate of the gene/CDS in genome.
        gene_end: 1-based end coordinate of the gene/CDS in genome.
        strand: '+' for positive strand, '-' for negative strand.
        aa_start: 1-based amino acid start position of the domain.
        aa_stop: 1-based amino acid end position of the domain.

    Returns:
        Tuple of (genome_start, genome_end) where genome_start <= genome_end.
        Returns (None, None) if strand is not '+' or '-'.
    """
    if strand == "+":
        g_start = gene_start + (aa_start - 1) * 3
        g_end = gene_start + (aa_stop * 3) - 1
        return g_start, g_end
    elif strand == "-":
        g_start = gene_end - (aa_stop * 3) + 1
        g_end = gene_end - (aa_start - 1) * 3
        return g_start, g_end
    else:
        return None, None
