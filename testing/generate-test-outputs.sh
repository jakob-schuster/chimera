TEST_INPUT=testing/second-run/CAG-B_S20_L001_R1_001/1000

./target/release/chimera \
    -i ${TEST_INPUT}/input.fastq \
    -r testing/second-run/reference.tsv \
    -o ${TEST_INPUT}/careful-output.tsv \
    --careful

./target/release/chimera \
    -i ${TEST_INPUT}/input.fastq \
    -r testing/second-run/reference.tsv \
    -o ${TEST_INPUT}/quick-output.tsv