+----------------+       +-------------------+       +-------------------+

| Sequencing |------>| Raw Data |------>| FASTQ File |
| Instrument | | (Reads) | | (Sequences + |
| | | | | Qualities) |
+----------------+       +-------------------+       +-------------------+
|
|
                         (Alignment to Reference)
|
                                   v
+----------------+       +-------------------+       +-------------------+

| Reference |------>| Aligned Data |------>| SAM/BAM/CRAM |
| Genome | | (Mapped Reads) | | File |
+----------------+       +-------------------+       +-------------------+
|
|
                         (Variant Calling)
|
                                   v
+----------------+       +-------------------+

| Variant Calls |------>| VCF File |
| (SNPs, etc.) | | (List of |
+----------------+ | Variants) |
                         +-------------------+
