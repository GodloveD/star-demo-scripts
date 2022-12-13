show_help() {
    echo "${0} < --all >"
}

if [ $# -gt 1 ]; then
    echo "too many arguments"
    show_help
    exit 1 
fi

ALL=0
if [ $# -eq 1 ]; then
    if [ "${1}x" != "--allx" ]; then 
        echo "unknown argument"
        show_help
        exit 1
    else 
        ALL=1
        echo "downloading all files"
    fi
fi

# download genome in fasta format
LINK_PREFIX="ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria//fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/"
LINK_LIST="Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.chromosome.Chromosome.fa.gz
Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz
Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.toplevel.fa.gz
Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna_rm.chromosome.Chromosome.fa.gz
Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna_rm.genome.fa.gz
Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna_rm.toplevel.fa.gz
Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna_sm.chromosome.Chromosome.fa.gz
Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna_sm.genome.fa.gz
Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna_sm.toplevel.fa.gz"

if [ $ALL -eq 1 ]; then
    for LINK in ${LINK_LIST}; do
        wget "${LINK_PREFIX}${LINK}"
    done
else
    LINK="Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz"
    wget "${LINK_PREFIX}${LINK}"
fi

# download annotations in gtf format
wget "ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria//gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf.gz"
