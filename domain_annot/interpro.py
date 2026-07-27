"""InterPro entry types, TSV parsing, and GO term mapping module."""

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, TextIO, Union


TYPE_PRIORITIES: Dict[str, int] = {
    "Family": 3,
    "Domain": 2,
    "Homologous_superfamily": 1,
}


def get_type_priority(entry_type: Optional[str]) -> int:
    """Return integer priority for InterPro entry type (Family=3, Domain=2, Homologous_superfamily=1, default=0)."""
    if not entry_type:
        return 0
    return TYPE_PRIORITIES.get(entry_type, 0)


@dataclass
class InterProEntryType:
    interpro_acc: str
    entry_type: str
    entry_name: str


@dataclass
class RawInterProHit:
    protein_id: str
    md5: str
    seq_len: int
    analysis: str
    signature_acc: str
    signature_desc: str
    start: int
    stop: int
    score: Optional[float]
    status: str
    date: str
    interpro_acc: Optional[str]
    interpro_desc: Optional[str]
    go_terms: Optional[str]
    pathways: Optional[str]


def parse_interpro_entry_list(
    file_path_or_handle: Union[str, Path, TextIO]
) -> Dict[str, InterProEntryType]:
    """
    Parse EBI interpro_entry_list.txt file.

    Columns: ENTRY_AC, ENTRY_TYPE, ENTRY_NAME
    Returns dict mapping interpro_acc -> InterProEntryType
    """
    entries: Dict[str, InterProEntryType] = {}

    if isinstance(file_path_or_handle, (str, Path)):
        handle = open(file_path_or_handle, "r")
        should_close = True
    else:
        handle = file_path_or_handle
        should_close = False

    try:
        reader = csv.reader(handle, delimiter="\t")
        first_line = True
        for row in reader:
            if not row or len(row) < 3:
                continue
            if first_line and ("ENTRY_AC" in row[0] or "interpro_acc" in row[0]):
                first_line = False
                continue
            first_line = False

            acc, etype, ename = row[0].strip(), row[1].strip(), row[2].strip()
            entries[acc] = InterProEntryType(
                interpro_acc=acc,
                entry_type=etype,
                entry_name=ename
            )
    finally:
        if should_close:
            handle.close()

    return entries


def parse_interpro_tsv(
    file_path_or_handle: Union[str, Path, TextIO]
) -> List[RawInterProHit]:
    """
    Parse InterProScan output TSV file (14 or 15 standard TSV columns).

    Columns:
    0: protein_id
    1: md5
    2: seq_len
    3: analysis
    4: signature_acc
    5: signature_desc
    6: start
    7: stop
    8: score
    9: status
    10: date
    11: interpro_acc
    12: interpro_desc
    13: go_terms
    14: pathways
    """
    hits: List[RawInterProHit] = []

    if isinstance(file_path_or_handle, (str, Path)):
        handle = open(file_path_or_handle, "r")
        should_close = True
    else:
        handle = file_path_or_handle
        should_close = False

    try:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or len(row) < 9:
                continue

            protein_id = row[0].strip()
            md5 = row[1].strip()
            seq_len = int(row[2])
            analysis = row[3].strip()
            signature_acc = row[4].strip()
            signature_desc = row[5].strip()
            start = int(row[6])
            stop = int(row[7])

            raw_score = row[8].strip()
            score: Optional[float] = None
            if raw_score and raw_score != "-" and raw_score.upper() != "NA":
                try:
                    score = float(raw_score)
                except ValueError:
                    score = None

            status = row[9].strip() if len(row) > 9 else "T"
            date = row[10].strip() if len(row) > 10 else ""
            
            interpro_acc = row[11].strip() if len(row) > 11 and row[11].strip() not in ("-", "", "NA") else None
            interpro_desc = row[12].strip() if len(row) > 12 and row[12].strip() not in ("-", "", "NA") else None
            go_terms = row[13].strip() if len(row) > 13 and row[13].strip() not in ("-", "", "NA") else None
            pathways = row[14].strip() if len(row) > 14 and row[14].strip() not in ("-", "", "NA") else None

            hits.append(RawInterProHit(
                protein_id=protein_id,
                md5=md5,
                seq_len=seq_len,
                analysis=analysis,
                signature_acc=signature_acc,
                signature_desc=signature_desc,
                start=start,
                stop=stop,
                score=score,
                status=status,
                date=date,
                interpro_acc=interpro_acc,
                interpro_desc=interpro_desc,
                go_terms=go_terms,
                pathways=pathways
            ))
    finally:
        if should_close:
            handle.close()

    return hits


def extract_go_mapping(hits: List[RawInterProHit]) -> Dict[str, str]:
    """
    Extract unique GO term IDs for each InterPro Accession,
    returning a dict mapping interpro_acc -> '|'-separated GO terms (e.g. 'GO:0003677|GO:0006355').
    """
    acc_to_go: Dict[str, Set[str]] = {}

    for hit in hits:
        if not hit.interpro_acc or not hit.go_terms:
            continue
        go_ids = re.findall(r"GO:\d+", hit.go_terms)
        if go_ids:
            if hit.interpro_acc not in acc_to_go:
                acc_to_go[hit.interpro_acc] = set()
            acc_to_go[hit.interpro_acc].update(go_ids)

    return {acc: "|".join(sorted(list(go_set))) for acc, go_set in acc_to_go.items()}


@dataclass
class InterProParentInfo:
    parent_acc: str
    parent_name: str


def parse_parent_child_tree(
    file_path_or_handle: Union[str, Path, TextIO]
) -> Dict[str, InterProParentInfo]:
    """
    Parse EBI ParentChildTreeFile.txt file.

    Returns dict mapping child interpro_acc -> InterProParentInfo(parent_acc, parent_name).
    """
    parent_map: Dict[str, InterProParentInfo] = {}

    if isinstance(file_path_or_handle, (str, Path)):
        handle = open(file_path_or_handle, "r")
        should_close = True
    else:
        handle = file_path_or_handle
        should_close = False

    try:
        stack: Dict[int, InterProParentInfo] = {}
        for line in handle:
            line_str = line.strip()
            if not line_str:
                continue

            lstripped = line_str.lstrip("-")
            dash_count = len(line_str) - len(lstripped)
            depth = dash_count // 2

            parts = lstripped.split("::")
            if not parts or not parts[0].startswith("IPR"):
                continue

            acc = parts[0].strip()
            name = parts[1].strip() if len(parts) > 1 else ""

            if depth > 0 and (depth - 1) in stack:
                parent_info = stack[depth - 1]
                parent_map[acc] = parent_info

            stack[depth] = InterProParentInfo(parent_acc=acc, parent_name=name)

    finally:
        if should_close:
            handle.close()

    return parent_map
