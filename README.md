# Differential Gene Expression Analysis.

[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) can be used on any *POSIX* compatible system (Linux, OS X, etc). It requires **Bash 3.2** (or later) and **Java 8** (or later, up to 11) to be installed.

1. Have atleast **Java 8** or later installed. Check if **Java** is installed using the following command:

```
$ java --version
openjdk 11.0.5 2019-10-15
OpenJDK Runtime Environment (build 11.0.5+10-post-Ubuntu-2ubuntu116.04)
OpenJDK 64-Bit Server VM (build 11.0.5+10-post-Ubuntu-2ubuntu116.04, mixed mode, sharing)
```

2. Download the the **nextflow** main executable file using the following command: 

```
$ wget -qO- https://get.nextflow.io | bash

      N E X T F L O W
      version 20.01.0 build 5264
      created 12-02-2020 10:14 UTC (10:14 GMT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io


Nextflow installation completed. Please note:
- the executable file `nextflow` has been created in the folder: /mnt/c/Users/noel.namai/personal/asimov
- you may complete the installation by moving it to a directory in your $PATH
```

3. Move the **nextflow** executable file to a directory accessible by the **$PATH** variable e.g.

```
$ mv nextflow /usr/local/bin/
```

```
$ nextflow -version

      N E X T F L O W
      version 20.01.0 build 5264
      created 12-02-2020 10:14 UTC (10:14 GMT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

```
$ tree
├── data
│   ├── ecoli_dh10b_annotation.fasta
│   ├── ecoli_dh10b_annotation.gff3
│   ├── ecoli_dh10b_ensembl.fasta
│   ├── ecoli_dh10b_ensembl.gff3
│   ├── ecoli_dh10b_ncbi.fasta
│   ├── ecoli_state1_rep1.fastq.gz
│   ├── ecoli_state1_rep1_2.fastq.gz
│   ├── ecoli_state1_rep2.fastq.gz
│   ├── ecoli_state1_rep2_2.fastq.gz
│   ├── ecoli_state2_rep1.fastq.gz
│   ├── ecoli_state2_rep1_2.fastq.gz
│   ├── ecoli_state2_rep2.fastq.gz
│   ├── ecoli_state2_rep2_2.fastq.gz
│   └── plasmid.fasta
├── dgexp.nf
├── docker
│   ├── Dockerfile
│   └── dgexp.R
├── nextflow.config
└── results
    ├── dexp-genes.tsv
    ├── genes-results.tsv
    ├── ma-plot.png
    ├── max-counts-plot.png
    └── min-counts-plot.png
```
