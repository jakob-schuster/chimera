TEST_INPUT=testing/second-run/CAG-B_S20_L001_R1_001/1000

./target/release/chimera \
    -r testing/second-run/reference.tsv \
    -i ${TEST_INPUT}/input.fastq \
    -c ${TEST_INPUT}/careful-chimera.fastq \
    -v ${TEST_INPUT}/careful-valid.fastq \
    -o ${TEST_INPUT}/careful-output.tsv \
    --careful

./target/release/chimera \
    -r testing/second-run/reference.tsv \
    -i ${TEST_INPUT}/input.fastq \
    -c ${TEST_INPUT}/quick-chimera.fastq \
    -v ${TEST_INPUT}/quick-valid.fastq \
    -o ${TEST_INPUT}/quick-output.tsv