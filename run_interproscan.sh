#!/usr/bin/env bash

# Required (placed in data_dir by the user; not downloaded by this script):
# - Mabs_proteins.nostop.fasta (reference strain protein sequences, no asterisks for stop codons)
# - Mabs.gff3 (reference strain GFF annotation file)


wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.77-108.0/interproscan-5.77-108.0-64-bit.tar.gz
tar -xzvf interproscan-5.77-108.0-64-bit.tar.gz
cd interproscan-5.77-108.0


./interproscan.sh \
    -f TSV \
    -i ${data_dir}/Mabs_proteins.nostop.fasta  \
    -iprlookup -goterms -pathways -f tsv \
    -d ${work_dir}/interproscan_results/ \
    -t p -T ${work_dir}/tmp/ -cpu 32


# list of domains:
wget ftp://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list \
    -O ${data_dir}/interpro_entry_list.txt
