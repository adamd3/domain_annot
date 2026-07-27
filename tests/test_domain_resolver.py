"""Unit tests for domain_resolver module."""

from domain_annot.gff import GeneFeature
from domain_annot.interpro import InterProEntryType, RawInterProHit
from domain_annot.domain_resolver import resolve_domains


def test_domain_resolution_tetr():
    entry_list = {
        "IPR001647": InterProEntryType("IPR001647", "Domain", "DNA-binding HTH domain, TetR-type"),
        "IPR009057": InterProEntryType("IPR009057", "Homologous_superfamily", "Homeodomain-like superfamily"),
        "IPR023772": InterProEntryType("IPR023772", "Conserved_site", "TetR conserved site"),
        "IPR036271": InterProEntryType("IPR036271", "Homologous_superfamily", "Tetracyclin repressor C-terminal"),
        "IPR050109": InterProEntryType("IPR050109", "Family", "TetR transcriptional regulator family"),
    }

    hits = [
        # Group 1: IPR009057 (Superfamily) 15-87
        RawInterProHit("MAB_RS20425", "md5", 219, "SUPERFAMILY", "SSF46689", "Homeo", 15, 87, 6.48e-19, "T", "", "IPR009057", "Homeodomain-like superfamily", None, None),
        # Group 2: IPR023772 (Conserved site) 39-70
        RawInterProHit("MAB_RS20425", "md5", 219, "ProSitePatterns", "PS01081", "TetR", 39, 70, None, "T", "", "IPR023772", "TetR conserved site", None, None),
        # Group 3: IPR050109 (Family) 11-211 (Full length parent -> should be filtered out)
        RawInterProHit("MAB_RS20425", "md5", 219, "PANTHER", "PTHR30055", "RUTR", 11, 211, 2.20e-18, "T", "", "IPR050109", "TetR transcriptional regulator family", None, None),
        # Group 4: IPR036271 (Superfamily) 102-183
        RawInterProHit("MAB_RS20425", "md5", 219, "SUPERFAMILY", "SSF48498", "TetC", 102, 183, 5.34e-07, "T", "", "IPR036271", "Tetracyclin repressor C-terminal", None, None),
        # Group 5: IPR001647 (Domain) 28-73 and 21-81 -> min start 21, max stop 81
        RawInterProHit("MAB_RS20425", "md5", 219, "Pfam", "PF00440", "TetR", 28, 73, 6.60e-11, "T", "", "IPR001647", "DNA-binding HTH domain, TetR-type", "GO:0003677(InterPro)", None),
        RawInterProHit("MAB_RS20425", "md5", 219, "ProSiteProfiles", "PS50977", "TetR", 21, 81, 22.45, "T", "", "IPR001647", "DNA-binding HTH domain, TetR-type", "GO:0003677(InterPro)", None),
    ]

    genes = {
        "MAB_RS20425": GeneFeature("MAB_RS20425", "chr1", 4072028, 4072684, "+")
    }

    resolved = resolve_domains(hits, entry_list, genes)

    # Should resolve 2 distinct domains:
    # 1. N-terminal HTH domain (IPR001647 or specificity upgraded)
    # 2. C-terminal domain (IPR036271)
    assert len(resolved) == 2
    accs = [r.interpro_acc for r in resolved]
    assert "IPR001647" in accs or "IPR009057" in accs
    assert "IPR036271" in accs

    # Check coordinate mapping on positive strand for resolved domain
    d1 = resolved[0]
    assert d1.strand == "+"
    assert d1.chrom == "chr1"
    assert d1.genome_start is not None
    assert d1.genome_end is not None
    assert d1.genome_start < d1.genome_end


def test_parent_domain_resolution():
    from domain_annot.interpro import InterProParentInfo, parse_parent_child_tree
    import io

    sample_tree = """IPR003593::AAA+ ATPase domain::
--IPR020591::Chromosomal replication control, initiator DnaA-like::
"""
    parent_map = parse_parent_child_tree(io.StringIO(sample_tree))
    assert "IPR020591" in parent_map
    assert parent_map["IPR020591"].parent_acc == "IPR003593"
    assert parent_map["IPR020591"].parent_name == "AAA+ ATPase domain"

    entry_list = {
        "IPR020591": InterProEntryType("IPR020591", "Family", "Chromosomal replication control, initiator DnaA-like")
    }
    hits = [
        RawInterProHit("prot1", "md5", 500, "PRINTS", "PR00051", "DnaA", 212, 491, 1e-50, "T", "", "IPR020591", "Chromosomal replication control, initiator DnaA-like", None, None)
    ]
    resolved = resolve_domains(hits, entry_list, parent_map=parent_map)
    assert len(resolved) == 1
    assert resolved[0].parent_acc == "IPR003593"
    assert resolved[0].parent_name == "AAA+ ATPase domain"

