"""Exporters for TSV, BED, and GFF3 domain annotation output files."""

import csv
from pathlib import Path
from typing import List, Union
from domain_annot.domain_resolver import ResolvedDomain


def write_tsv(domains: List[ResolvedDomain], out_path: Union[str, Path]) -> Path:
    """
    Write resolved domain annotations to a TSV file.

    Columns:
    protein_id, aa_start, aa_stop, strand, genome_start, genome_end, interpro_acc, interpro_desc, entry_type, go_terms
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    headers = [
        "protein_id", "aa_start", "aa_stop", "strand",
        "genome_start", "genome_end", "interpro_acc",
        "interpro_desc", "entry_type", "parent_acc",
        "parent_name", "go_terms"
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)

        for d in domains:
            writer.writerow([
                d.protein_id,
                d.aa_start,
                d.aa_stop,
                d.strand or "",
                d.genome_start if d.genome_start is not None else "",
                d.genome_end if d.genome_end is not None else "",
                d.interpro_acc,
                d.interpro_desc,
                d.entry_type or "",
                d.parent_acc or "",
                d.parent_name or "",
                d.go_terms or ""
            ])

    return out_path


def write_bed(domains: List[ResolvedDomain], out_path: Union[str, Path]) -> Path:
    """
    Write resolved domain annotations to 6-column BED format.

    Columns:
    chrom, chromStart (0-based), chromEnd (1-based), name, score, strand
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as f:
        for d in domains:
            if not d.chrom or d.genome_start is None or d.genome_end is None:
                continue
            # BED coordinates: chromStart is 0-based, chromEnd is 1-based exclusive
            chrom_start = d.genome_start - 1
            chrom_end = d.genome_end
            name = f"{d.protein_id}:{d.interpro_acc}"
            strand = d.strand or "+"
            f.write(f"{d.chrom}\t{chrom_start}\t{chrom_end}\t{name}\t0\t{strand}\n")

    return out_path


def write_gff3(domains: List[ResolvedDomain], out_path: Union[str, Path]) -> Path:
    """
    Write resolved domain annotations to GFF3 format.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as f:
        f.write("##gff-version 3\n")
        for i, d in enumerate(domains, 1):
            if not d.chrom or d.genome_start is None or d.genome_end is None:
                continue
            chrom = d.chrom
            source = "DomainAnnot"
            ftype = "protein_domain"
            start = d.genome_start
            end = d.genome_end
            score = "."
            strand = d.strand or "+"
            phase = "."
            
            # Attributes
            attr_parts = [
                f"ID=domain-{d.protein_id}-{i}",
                f"Name={d.interpro_acc}",
                f"Target={d.protein_id}",
                f"Dbxref=InterPro:{d.interpro_acc}",
                f"Note={d.interpro_desc}"
            ]
            if d.entry_type:
                attr_parts.append(f"entry_type={d.entry_type}")
            if d.parent_acc:
                attr_parts.append(f"parent_acc={d.parent_acc}")
            if d.parent_name:
                attr_parts.append(f"parent_name={d.parent_name}")
            if d.go_terms:
                attr_parts.append(f"Ontology_term={d.go_terms}")

            attrs = ";".join(attr_parts)
            f.write(f"{chrom}\t{source}\t{ftype}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attrs}\n")

    return out_path
