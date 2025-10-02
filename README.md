# fastq_to_fasta

This is a small experimental script I wrote to test the Galaxy API. It uploads one or two FASTQ files to a Galaxy instance, tries to run a de novo RNA assembler (like rnaSPAdes or Trinity), and downloads the resulting FASTA file. If no assembler is available, it falls back to converting FASTQ to FASTA. This was just one of my experiments exploring automation with Galaxy workflows.
