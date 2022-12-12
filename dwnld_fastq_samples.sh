show_help() {
    echo "${0} < integer file number >"
}

if [ $# -gt 1 ]; then
    echo "too many arguments"
    show_help
    exit 1
fi

FASTQ_FILES="SRX10348166
SRX10348167
SRX10348168
SRX10348169
SRX10348170
SRX10348171
SRX10348172
SRX10348173
SRX10348174
SRX10348175
SRX10348176
SRX10348177
SRX10348178
SRX10348179
SRX10348180
SRX10348181
"

if [ $# -eq 1 ]; then
    if [[ "$1" =~ ^[0-9]+$ ]]; then
        FASTQ_FILES=$(echo "${FASTQ_FILES}" | head -n $1 | tail -n 1)
    else
        echo "unknown argument $1"
        show_help
        exit 1
    fi
else
    echo "downloading all files"
fi

if [ "${FASTQ_FILES}X" == "X" ]; then 
    echo "integer out of bounds"
    exit 1
fi

for FILE in $FASTQ_FILES; do
    wget -O "${FILE}.fq.gz" "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=${FILE}"
done



