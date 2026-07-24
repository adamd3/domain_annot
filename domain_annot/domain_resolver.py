"""Domain resolution engine: boundary aggregation, specificity upgrade, overlap reduction, and coordinate mapping."""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional

from domain_annot.coordinate_mapper import map_aa_to_genome_coords
from domain_annot.gff import GeneFeature
from domain_annot.interpro import (
    InterProEntryType,
    RawInterProHit,
    extract_go_mapping,
    get_type_priority,
)


@dataclass
class ResolvedDomain:
    protein_id: str
    aa_start: int
    aa_stop: int
    interpro_acc: str
    interpro_desc: str
    entry_type: Optional[str]
    best_evalue: Optional[float]
    go_terms: Optional[str] = None
    strand: Optional[str] = None
    genome_start: Optional[int] = None
    genome_end: Optional[int] = None
    chrom: Optional[str] = None


@dataclass
class _CandidateDomain:
    protein_id: str
    start: int
    stop: int
    interpro_acc: str
    interpro_desc: str
    entry_type: Optional[str]
    type_priority: int
    best_evalue: Optional[float]

    @property
    def width(self) -> int:
        return self.stop - self.start + 1


def resolve_domains(
    hits: List[RawInterProHit],
    entry_list: Dict[str, InterProEntryType],
    genes: Optional[Dict[str, GeneFeature]] = None
) -> List[ResolvedDomain]:
    """
    Resolve raw InterProScan hits into non-redundant domain annotations mapped to genomic coordinates.

    Args:
        hits: List of RawInterProHit objects from InterProScan output.
        entry_list: Dict mapping interpro_acc -> InterProEntryType from EBI entry list.
        genes: Optional dict mapping protein_id -> GeneFeature for genomic coordinate translation.

    Returns:
        List of ResolvedDomain objects sorted by protein_id and aa_start.
    """
    # 1. Extract GO mapping lookup per interpro_acc
    go_map = extract_go_mapping(hits)

    # 2. Group raw hits by (protein_id, interpro_acc) to calculate outermost start/stop & best E-value
    grouped_hits: Dict[tuple[str, str], List[RawInterProHit]] = {}
    for hit in hits:
        if not hit.interpro_acc:
            continue
        key = (hit.protein_id, hit.interpro_acc)
        if key not in grouped_hits:
            grouped_hits[key] = []
        grouped_hits[key].append(hit)

    # Build initial candidate boundaries
    protein_candidates: Dict[str, List[_CandidateDomain]] = {}
    for (protein_id, interpro_acc), group in grouped_hits.items():
        min_start = min(h.start for h in group)
        max_stop = max(h.stop for h in group)
        
        # Determine best E-value
        valid_scores = [h.score for h in group if h.score is not None]
        best_evalue = min(valid_scores) if valid_scores else None

        # Determine interpro description
        desc_list = [h.interpro_desc for h in group if h.interpro_desc]
        interpro_desc = desc_list[0] if desc_list else ""

        # Lookup entry type and priority
        entry_meta = entry_list.get(interpro_acc)
        entry_type = entry_meta.entry_type if entry_meta else None
        type_priority = get_type_priority(entry_type)

        cand = _CandidateDomain(
            protein_id=protein_id,
            start=min_start,
            stop=max_stop,
            interpro_acc=interpro_acc,
            interpro_desc=interpro_desc,
            entry_type=entry_type,
            type_priority=type_priority,
            best_evalue=best_evalue
        )

        if protein_id not in protein_candidates:
            protein_candidates[protein_id] = []
        protein_candidates[protein_id].append(cand)

    # 3. Process candidates per protein (filter parent structural entries, upgrade specificity, drop >50% overlaps)
    final_domains: List[ResolvedDomain] = []

    for protein_id, candidates in protein_candidates.items():
        if not candidates:
            continue

        widths = [c.width for c in candidates]
        max_protein_width = max(widths)
        num_candidates = len(candidates)

        # Filter out full-length parent entries (width == max_protein_width and width > 150 and num_candidates > 1)
        filtered = [
            c for c in candidates
            if not (c.width == max_protein_width and c.width > 150 and num_candidates > 1)
        ]
        if not filtered:
            filtered = candidates

        # Sort by width descending
        sorted_cands = sorted(filtered, key=lambda c: c.width, reverse=True)

        claimed_segments: List[tuple[int, int]] = []
        
        for cand in sorted_cands:
            # Check for nested child hits to perform Specificity Upgrade
            nested = [
                c for c in sorted_cands
                if c.start >= cand.start and c.stop <= cand.stop and c.interpro_acc != cand.interpro_acc
            ]

            if nested:
                # Sort nested by desc(type_priority), then best_evalue
                def sort_key(c: _CandidateDomain):
                    eval_val = c.best_evalue if c.best_evalue is not None else float("inf")
                    return (-c.type_priority, eval_val)

                top_child = sorted(nested, key=sort_key)[0]
                if top_child.type_priority > cand.type_priority:
                    # Inherit child identity
                    cand.interpro_acc = top_child.interpro_acc
                    cand.interpro_desc = top_child.interpro_desc
                    cand.entry_type = top_child.entry_type
                    cand.type_priority = top_child.type_priority

            # Overlap Reduction Check against claimed segments
            is_redundant = False
            cand_width = cand.width

            for seg_start, seg_stop in claimed_segments:
                overlap_start = max(cand.start, seg_start)
                overlap_stop = min(cand.stop, seg_stop)

                if overlap_start <= overlap_stop:
                    overlap_width = overlap_stop - overlap_start + 1
                    if (overlap_width / cand_width) > 0.5:
                        is_redundant = True
                        break

            if not is_redundant:
                claimed_segments.append((cand.start, cand.stop))

                # Lookup GO terms
                go_terms = go_map.get(cand.interpro_acc)

                # Translate genomic coordinates if gene mapping provided
                strand: Optional[str] = None
                genome_start: Optional[int] = None
                genome_end: Optional[int] = None
                chrom: Optional[str] = None

                if genes and protein_id in genes:
                    gf = genes[protein_id]
                    strand = gf.strand
                    chrom = gf.chrom
                    genome_start, genome_end = map_aa_to_genome_coords(
                        gene_start=gf.start,
                        gene_end=gf.end,
                        strand=gf.strand,
                        aa_start=cand.start,
                        aa_stop=cand.stop
                    )

                final_domains.append(ResolvedDomain(
                    protein_id=protein_id,
                    aa_start=cand.start,
                    aa_stop=cand.stop,
                    interpro_acc=cand.interpro_acc,
                    interpro_desc=cand.interpro_desc,
                    entry_type=cand.entry_type,
                    best_evalue=cand.best_evalue,
                    go_terms=go_terms,
                    strand=strand,
                    genome_start=genome_start,
                    genome_end=genome_end,
                    chrom=chrom
                ))

    # Sort final domains by protein_id and aa_start
    final_domains.sort(key=lambda d: (d.protein_id, d.aa_start))
    return final_domains
