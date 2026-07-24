# eval "$(/home/ad866/miniforge3/bin/conda shell.bash hook)"
# mamba activate new_bact_seq



library(tidyverse)
library(vroom)
library(RColorBrewer)
library(IRanges)


work_dir <- "/scratch/ad866/Mabs_proteomics"

fname <- "parse_interpro_scan.RData"
load(file = file.path(work_dir, fname))



########################
## Colour palettes
########################
large_disc_pal <- brewer.pal.info[brewer.pal.info$category == "qual", ]
colpal_large <- unlist(
    mapply(brewer.pal, large_disc_pal$maxcolors, rownames(large_disc_pal))
)
colpal_large[c(5:8)] <- colpal_large[c(70:73)] ## replace to avoid colour clashes


gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
ggColsDefault <- (gg_color_hue(4))


## -----------------------------------------------------------------------------
## Read data
## -----------------------------------------------------------------------------

ref_gff <- read_tsv(
    file.path(work_dir, "Mabs.gff3"),
    comment = "#", col_names = FALSE
) %>%
    filter(X3 == "gene") %>% # Only look at genes
    select(chrom = X1, start = X4, end = X5, strand = X7, info = X9) %>%
    mutate(gene_name = gsub(
        "gene-", "",
        str_extract(info, "(?<=Name=)[^;]+|(?<=ID=)[^;]+")
    )) %>%
    select(gene_name, start, end, strand)

ref_gff %>% head()
# # A tibble: 6 × 4
#   gene_name   start   end strand
#   <chr>       <dbl> <dbl> <chr>
# 1 MAB_RS00150     1  1476 +
# 2 MAB_RS00155  1553  1747 -
# 3 MAB_RS00160  2169  3368 +
# 4 MAB_RS00165  3419  4312 +
# 5 MAB_RS00170  4332  5480 +
# 6 MAB_RS00175  5473  6024 +


# Load the universal EBI entry classifications for InterPro Accessions
ipr_types <- vroom(
    file.path(
        work_dir,
        "interpro_entry_list.txt"
    ),
    delim = "\t",
    skip = 1,
    col_names = c("interpro_acc", "entry_type", "entry_name"),
    show_col_types = FALSE
) %>%
    select(interpro_acc, entry_type)

# count entries per type:
ipr_types %>%
    group_by(entry_type) %>%
    summarize(n = n()) %>%
    data.frame()
#               entry_type     n
# 1            Active_site   133
# 2           Binding_site    82
# 3         Conserved_site   775
# 4                 Domain 21357
# 5                 Family 27926
# 6 Homologous_superfamily  3510
# 7                    PTM    17
# 8                 Repeat   390


# Load the raw InterProScan results for our M. abscessus reference proteome
ipr_columns <- c(
    "protein_id", "md5", "seq_len", "analysis", "signature_acc",
    "signature_desc", "start", "stop", "score", "status",
    "date", "interpro_acc", "interpro_desc", "go_terms", "pathways"
)

ipr_raw <- vroom(
    file.path(
        work_dir,
        "Mabs_proteins.nostop.fasta.tsv"
    ),
    delim = "\t",
    col_names = ipr_columns,
    na = c("-", "", "NA"), # Automatically turns hyphens into real NAs
    show_col_types = FALSE
)

ipr_raw %>%
    select(-pathways) %>%
    filter(protein_id == "MAB_RS20425") %>%
    data.frame()

#     protein_id                              md5 seq_len        analysis
# 1  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219     SUPERFAMILY
# 2  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219 ProSitePatterns
# 3  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219         PANTHER
# 4  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219     SUPERFAMILY
# 5  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219            Pfam
# 6  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219          PRINTS
# 7  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219          PRINTS
# 8  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219          Gene3D
# 9  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219      MobiDBLite
# 10 MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219 ProSiteProfiles
#        signature_acc                                  signature_desc start stop
# 1           SSF46689                                Homeodomain-like    15   87
# 2            PS01081                 TetR-type HTH domain signature.    39   70
# 3          PTHR30055         HTH-TYPE TRANSCRIPTIONAL REGULATOR RUTR    11  211
# 4           SSF48498   Tetracyclin repressor-like, C-terminal domain   102  183
# 5            PF00440      Bacterial regulatory proteins, tetR family    28   73
# 6            PR00455 TetR bacterial regulatory protein HTH signature    48   71
# 7            PR00455 TetR bacterial regulatory protein HTH signature    27   40
# 8  G3DSA:1.10.357.10                Tetracycline Repressor, domain 2    15  211
# 9        mobidb-lite                   consensus disorder prediction     1   22
# 10           PS50977                   TetR-type HTH domain profile.    21   81
#          score status       date interpro_acc
# 1  6.48000e-19   TRUE 17-06-2026    IPR009057
# 2           NA   TRUE 17-06-2026    IPR023772
# 3  2.20000e-18   TRUE 17-06-2026    IPR050109
# 4  5.34000e-07   TRUE 17-06-2026    IPR036271
# 5  6.60000e-11   TRUE 17-06-2026    IPR001647
# 6  9.50000e-08   TRUE 17-06-2026    IPR001647
# 7  9.50000e-08   TRUE 17-06-2026    IPR001647
# 8  7.70000e-35   TRUE 17-06-2026         <NA>
# 9           NA   TRUE 17-06-2026         <NA>
# 10 2.24556e+01   TRUE 17-06-2026    IPR001647
#                                                interpro_desc
# 1                                Homedomain-like superfamily
# 2          DNA-binding HTH domain, TetR-type, conserved site
# 3              HTH-type, TetR-like transcriptional regulator
# 4  Tetracyclin repressor-like, C-terminal domain superfamily
# 5                          DNA-binding HTH domain, TetR-type
# 6                          DNA-binding HTH domain, TetR-type
# 7                          DNA-binding HTH domain, TetR-type
# 8                                                       <NA>
# 9                                                       <NA>
# 10                         DNA-binding HTH domain, TetR-type
#                                                                                                 go_terms
# 1                                                                                                   <NA>
# 2                                                                                                   <NA>
# 3  GO:0000976(PANTHER)|GO:0003700(InterPro)|GO:0003700(PANTHER)|GO:0006355(PANTHER)|GO:0006355(InterPro)
# 4                                                                                                   <NA>
# 5                                                                                   GO:0003677(InterPro)
# 6                                                                                   GO:0003677(InterPro)
# 7                                                                                   GO:0003677(InterPro)
# 8                                                                                                   <NA>
# 9                                                                                                   <NA>
# 10                                                                                  GO:0003677(InterPro)


# Create a clean lookup dictionary of InterPro Accessions -> GO Terms
interpro_go_map <- ipr_raw %>%
    filter(!is.na(interpro_acc) & !is.na(go_terms)) %>%
    distinct(interpro_acc, go_terms) %>%
    mutate(go_id = str_extract_all(go_terms, "GO:\\d+")) %>%
    unnest(go_id) %>%
    distinct(interpro_acc, go_id) %>%
    # Collapse them into a clean, comma-separated string per InterPro ID
    group_by(interpro_acc) %>%
    summarize(go_terms = paste(go_id, collapse = "|"), .groups = "drop")


## Get a clean, non-overlapping list of domains present in each protein
# The Strategy
# - Group by Protein & InterPro ID: Instead of collapsing everything
#   globally, we group overlapping coordinates that share the same
#   overarching biological concept (interpro_acc). For example,
#   rows 5 and 6 share IPR001647, so their coordinates (28–73 and 48–71)
#   should be merged into a single outer boundary (28–73).

# - Resolve Multi-Domain Architectures: For a TetR regulator, you
#   expect an N-terminal DNA-binding domain and a C-terminal
#   dimerization/ligand-binding domain.

# - Hierarchy Flattening: If two distinct InterPro IDs still overlap,
#   we prioritize them by selecting the entry with the best E-value or
#   the more granular description.

# Our data shows this perfectly:

# Domain 1 (N-terminus): IPR001647 / IPR023772 / IPR009057 cluster
# roughly between residues 15 and 87.
# Domain 2 (C-terminus): IPR036271 clusters between residues 102 and 183.


ipr_raw %>%
    filter(protein_id == "MAB_RS20425") %>%
    select(-go_terms, -pathways) %>%
    data.frame()
#     protein_id                              md5 seq_len        analysis
# 1  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219     SUPERFAMILY
# 2  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219 ProSitePatterns
# 3  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219         PANTHER
# 4  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219     SUPERFAMILY
# 5  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219            Pfam
# 6  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219          PRINTS
# 7  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219          PRINTS
# 8  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219          Gene3D
# 9  MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219      MobiDBLite
# 10 MAB_RS20425 cd391e55172e910bb6d3e816c5fe3a7b     219 ProSiteProfiles
#        signature_acc                                  signature_desc start stop
# 1           SSF46689                                Homeodomain-like    15   87
# 2            PS01081                 TetR-type HTH domain signature.    39   70
# 3          PTHR30055         HTH-TYPE TRANSCRIPTIONAL REGULATOR RUTR    11  211
# 4           SSF48498   Tetracyclin repressor-like, C-terminal domain   102  183
# 5            PF00440      Bacterial regulatory proteins, tetR family    28   73
# 6            PR00455 TetR bacterial regulatory protein HTH signature    48   71
# 7            PR00455 TetR bacterial regulatory protein HTH signature    27   40
# 8  G3DSA:1.10.357.10                Tetracycline Repressor, domain 2    15  211
# 9        mobidb-lite                   consensus disorder prediction     1   22
# 10           PS50977                   TetR-type HTH domain profile.    21   81
#          score status       date interpro_acc
# 1  6.48000e-19   TRUE 17-06-2026    IPR009057
# 2           NA   TRUE 17-06-2026    IPR023772
# 3  2.20000e-18   TRUE 17-06-2026    IPR050109
# 4  5.34000e-07   TRUE 17-06-2026    IPR036271
# 5  6.60000e-11   TRUE 17-06-2026    IPR001647
# 6  9.50000e-08   TRUE 17-06-2026    IPR001647
# 7  9.50000e-08   TRUE 17-06-2026    IPR001647
# 8  7.70000e-35   TRUE 17-06-2026         <NA>
# 9           NA   TRUE 17-06-2026         <NA>
# 10 2.24556e+01   TRUE 17-06-2026    IPR001647
#                                                interpro_desc
# 1                                Homedomain-like superfamily
# 2          DNA-binding HTH domain, TetR-type, conserved site
# 3              HTH-type, TetR-like transcriptional regulator
# 4  Tetracyclin repressor-like, C-terminal domain superfamily
# 5                          DNA-binding HTH domain, TetR-type
# 6                          DNA-binding HTH domain, TetR-type
# 7                          DNA-binding HTH domain, TetR-type
# 8                                                       <NA>
# 9                                                       <NA>
# 10                         DNA-binding HTH domain, TetR-type

# For each unique combination of protein_id + interpro_acc, calculate the
# outermost start and stop coordinates, and the best E-value
resolved_boundaries <- ipr_raw %>%
    filter(!is.na(interpro_acc)) %>%
    mutate(score = as.numeric(score)) %>%
    group_by(protein_id, interpro_acc, interpro_desc) %>%
    summarize(
        start = min(start),
        stop = max(stop),
        best_evalue = if (all(is.na(score))) NA_real_ else min(score, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    # Append the official EBI entry type to every row
    left_join(ipr_types, by = "interpro_acc") %>%
    # Assign numerical priority
    mutate(type_priority = case_when(
        entry_type == "Family" ~ 3,
        entry_type == "Domain" ~ 2,
        entry_type == "Homologous_superfamily" ~ 1,
        TRUE ~ 0
    ))

resolved_boundaries %>%
    filter(protein_id == "MAB_RS20425") %>%
    data.frame()
#    protein_id interpro_acc
# 1 MAB_RS20425    IPR001647
# 2 MAB_RS20425    IPR009057
# 3 MAB_RS20425    IPR023772
# 4 MAB_RS20425    IPR036271
# 5 MAB_RS20425    IPR050109
#                                               interpro_desc start stop
# 1                         DNA-binding HTH domain, TetR-type    21   81
# 2                               Homedomain-like superfamily    15   87
# 3         DNA-binding HTH domain, TetR-type, conserved site    39   70
# 4 Tetracyclin repressor-like, C-terminal domain superfamily   102  183
# 5             HTH-type, TetR-like transcriptional regulator    11  211
#   best_evalue             entry_type type_priority
# 1    6.60e-11                 Domain             2
# 2    6.48e-19 Homologous_superfamily             1
# 3          NA         Conserved_site             0
# 4    5.34e-07 Homologous_superfamily             1
# 5    2.20e-18                 Family             3


# Now we have a clean set of non-redundant domain boundaries per protein, but
# some of them still overlap.
# To resolve this, if one domain is fully nested within another, we check their
# type priorities. If the nested domain has a higher priority
# (e.g. Family > Domain > Homologous_superfamily), we inherit the identity of
# the nested domain for the parent domain. This allows us to capture
# multi-domain architectures without losing important annotations.

final_domains <- resolved_boundaries %>%
    mutate(best_evalue = ifelse(
        is.infinite(best_evalue), NA_real_, best_evalue
    )) %>%
    group_by(protein_id) %>%
    mutate(
        width = stop - start + 1,
        max_protein_width = max(width)
    ) %>%
    # Filter out full-length structural parent entries, but keep for
    # small proteins < 150aa [which likely only have a single domain, even
    # if there are multiple interpro hits] or single-domain proteins:
    filter(!(width == max_protein_width & width > 150 & n() > 1)) %>%
    arrange(desc(width)) %>%
    do({
        df <- .
        if (nrow(df) <= 1) {
            df
        } else {
            to_drop <- c()
            claimed_segments <- list()

            for (i in 1:nrow(df)) {
                current_start <- df$start[i]
                current_stop <- df$stop[i]
                current_width <- df$width[i]

                # Identify nested child hits inside this footprint
                nested_hits <- df %>%
                    filter(start >= current_start & stop <= current_stop &
                        interpro_acc != df$interpro_acc[i])

                if (nrow(nested_hits) > 0) {
                    # AUTOMATED UPGRADE: Sort nested hits by
                    # their official EBI Type Priority
                    specific_child <- nested_hits %>%
                        arrange(desc(type_priority), best_evalue) %>%
                        dplyr::slice(1)

                    # If the nested child is more specific than the
                    # parent, inherit its identity
                    if (nrow(specific_child) == 1 &&
                        specific_child$type_priority > df$type_priority[i]) {
                        df$interpro_acc[i] <- specific_child$interpro_acc
                        df$interpro_desc[i] <- specific_child$interpro_desc
                        df$entry_type[i] <- specific_child$entry_type
                        df$type_priority[i] <- specific_child$type_priority
                    }
                }

                # Standard Overlap Reduction Check
                is_redundant <- FALSE
                if (length(claimed_segments) > 0) {
                    for (seg in claimed_segments) {
                        overlap_start <- max(current_start, seg$start)
                        overlap_stop <- min(current_stop, seg$stop)

                        if (overlap_start <= overlap_stop) {
                            overlap_width <- overlap_stop - overlap_start + 1
                            if ((overlap_width / current_width) > 0.5) {
                                is_redundant <- TRUE
                                break
                            }
                        }
                    }
                }

                if (is_redundant) {
                    to_drop <- c(to_drop, i)
                } else {
                    claimed_segments[[length(
                        claimed_segments
                    ) + 1]] <- list(
                        start = current_start, stop = current_stop
                    )
                }
            }

            if (length(to_drop) > 0) df <- df[-to_drop, ]
            df
        }
    }) %>%
    ungroup() %>%
    mutate(width = stop - start + 1) %>%
    arrange(protein_id, start) %>%
    select(protein_id, start, stop, interpro_acc, interpro_desc, entry_type)


final_domains %>%
    filter(protein_id == "MAB_RS20425") %>%
    data.frame()
#    protein_id start stop interpro_acc
# 1 MAB_RS20425    15   87    IPR001647
# 2 MAB_RS20425   102  183    IPR036271
#                                               interpro_desc
# 1                         DNA-binding HTH domain, TetR-type
# 2 Tetracyclin repressor-like, C-terminal domain superfamily
#               entry_type
# 1                 Domain
# 2 Homologous_superfamily

# count unique proteins with at least one InterPro domain:
final_domains %>%
    pull(protein_id) %>%
    unique() %>%
    length()
# [1] 4275

# add the GO terms where available:
final_domains <- final_domains %>%
    left_join(interpro_go_map, by = "interpro_acc")

# count go terms assigned per entry_type:
final_domains %>%
    filter(!is.na(go_terms)) %>%
    group_by(entry_type) %>%
    summarize(n_with_go = n(), .groups = "drop") %>%
    arrange(desc(n_with_go)) %>%
    data.frame()
#               entry_type n_with_go
# 1                 Family      1208
# 2                 Domain      1153
# 3 Homologous_superfamily        94
# 4         Conserved_site        17
# 5            Active_site         2
# 6                 Repeat         1
# 7                   <NA>         1

# how many proteins have at least 1 GO term assigned via InterPro?
final_domains %>%
    filter(!is.na(go_terms)) %>%
    pull(protein_id) %>%
    unique() %>%
    length()
# [1] 2186

# Add genome nucleotide coordinates for each domain
final_domains <- final_domains %>%
    left_join(
        ref_gff %>% select(gene_name, gene_start = start, gene_end = end, strand),
        by = c("protein_id" = "gene_name")
    ) %>%
    # Perform strand-aware translation from AA coordinates to nt coordinates
    mutate(
        genome_start = case_when(
            strand == "+" ~ gene_start + (start - 1) * 3,
            strand == "-" ~ gene_end - (stop * 3) + 1,
            TRUE ~ NA_real_
        ),
        genome_end = case_when(
            strand == "+" ~ gene_start + (stop * 3) - 1,
            strand == "-" ~ gene_end - (start - 1) * 3,
            TRUE ~ NA_real_
        )
    ) %>%
    # Clean up and reorder columns for downstream usage (e.g., making a BED file)
    select(
        protein_id,
        aa_start = start,
        aa_stop = stop,
        strand,
        genome_start,
        genome_end,
        interpro_acc,
        interpro_desc,
        entry_type
    )

final_domains %>%
    filter(protein_id == "MAB_RS20425") %>%
    data.frame()
#    protein_id aa_start aa_stop strand genome_start genome_end interpro_acc
# 1 MAB_RS20425       15      87      +      4072070    4072288    IPR001647
# 2 MAB_RS20425      102     183      +      4072331    4072576    IPR036271
#                                               interpro_desc
# 1                         DNA-binding HTH domain, TetR-type
# 2 Tetracyclin repressor-like, C-terminal domain superfamily
#               entry_type
# 1                 Domain
# 2 Homologous_superfamily


final_domains %>%
    filter(protein_id == "MAB_RS00215") %>%
    data.frame()
#    protein_id aa_start aa_stop strand genome_start genome_end interpro_acc
# 1 MAB_RS00215        8     266      -        13119      13895    IPR001447
#                   interpro_desc entry_type
# 1 Arylamine N-acetyltransferase     Family
# ^ NOTE that for -ve strand genes, end coordinate is still greater than start


## -----------------------------------------------------------------------------
## Save workspace
## -----------------------------------------------------------------------------
save.image(file = file.path(work_dir, fname))
