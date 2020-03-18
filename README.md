# Differential Gene Expression Analysis

The final workflow, implemented using [Nextflow](https://www.nextflow.io/), can be found here: [dgexp.nf](https://github.com/noelnamai/dgexp/blob/master/dgexp.nf). The accompanying configuration file is available here: [nexflow.config](https://github.com/noelnamai/dgexp/blob/master/nextflow.config). 

The final table with differentially expressed genes is available here: [dexp-genes.tsv](https://github.com/noelnamai/dgexp/blob/master/results/dexp-genes.tsv). I used an *adjusted p value cutoff* of 0.1 and a *log2 foldchange* threshold for constructing **Wald tests of significance** of 1.5.

## Flow diagram

![Alt text](./data/dgexp.png)

## Software requirements

The **Dockerfile** used to generate the Docker container used in the workflow is available here: [Dockerfile](https://github.com/noelnamai/dgexp/blob/master/docker/Dockerfile).

Most of the tools used are standard off the shelf tools. However, **DESeq2** has been wrapped into an Rscript available here: [dgexp.R](https://github.com/noelnamai/dgexp/blob/master/docker/dgexp.R).

[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) can be used on any *POSIX* compatible system (Linux, OS X, etc). It requires **Bash 3.2** (or later) and **Java 8** (or later, up to 11) to be installed. It is important to run the workflow on a **Linux** based system with atleast **8 CPUs** and **30 GB of RAM**. This workflow was run and tested on an **m4.2xlarge Amazon Linux Instance**. 

## Work environment requirements

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

4. Install **Docker**. A good walk through is available here: [How To Install and Use Docker on Ubuntu 18.04](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04).

```
$ docker --version
Docker version 18.09.6, build 481bc77
```

5. Clone the github repository to your environment:

```
$ git clone https://github.com/noelnamai/dgexp.git
```

6. Set up the *data* directory to have the structure shown below. Both the **ecoli_dh10b_ensembl.fasta** and **ecoli_dh10b_ensembl.gff3** were downloaded from [Ensembl](http://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_dh10b/Info/Index). The provide **fasta** and **gff3** were not yielding any expression results for all the samples.

```
$ tree dgexp/

├── data
│   ├── ecoli_dh10b_ensembl.fasta
│   ├── ecoli_dh10b_ensembl.gff3
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

7. To run the workflow, use the following command from the **dgexp** directory:

```
$ nextflow run dgexp.nf

N E X T F L O W  ~  version 20.01.0
Launching `dgexp.nf` [spontaneous_kilby] - revision: 8693a87a25

D I F F E R E N T I A L  G E N E  E X P R E S I O N  A N A L Y S I S
====================================================================

Samtools   : 1.9
Gffread    : 0.11.8
HISAT2     : 2.1.0
HTSeq      : 0.11.1
DESeq2     : 1.26.0
Start time : 2020-03-18T15:17:01.419777Z

executor >  local (27)
[78/3e84e0] process > trim_reads             [100%] 4 of 4 ✔
[a7/89afb7] process > check_read_quality     [100%] 4 of 4 ✔
[68/53cdf5] process > convert_gff3_to_gtf    [100%] 1 of 1 ✔
[7b/0ad38b] process > extract_exons_and_ss   [100%] 1 of 1 ✔
[fb/fef2c1] process > build_genome_index     [100%] 1 of 1 ✔
[88/366215] process > map_reads_to_reference [100%] 4 of 4 ✔
[7f/cc3670] process > convert_sam_to_bam     [100%] 4 of 4 ✔
[65/d29197] process > sort_bam_file          [100%] 4 of 4 ✔
[69/959a52] process > generate_raw_counts    [100%] 4 of 4 ✔
[90/3dca4e] process > detect_dexp_genes      [100%] 1 of 1 ✔

Completed at: 18-Mar-2020 15:46:30
Duration    : 29m 29s
CPU hours   : 2.5
Succeeded   : 28
```

## Results

The following files are the final outputs from the *detect_dexp_genes* process which runs **DESeq2**:

1. **dexp-genes.tsv**: the table with the differentially expresed genes with *padj < 0.1* and *log2FoldChange > log2(1.5)*. The full table is available here [dexp-genes.tsv](https://github.com/noelnamai/dgexp/blob/master/results/dexp-genes.tsv).

|genes        |   baseMean| log2FoldChange|     lfcSE|     stat|    pvalue|      padj|
|:------------|----------:|--------------:|---------:|--------:|---------:|---------:|
|ECDH10B_0130 | 1220.43800|       2.525123| 0.4439391| 4.370330| 0.0000124| 0.0005456|
|ECDH10B_0664 |   69.08905|       5.171019| 1.1160349| 4.109241| 0.0000397| 0.0015977|
|ECDH10B_0739 | 1916.88143|       1.734633| 0.4012711| 2.865072| 0.0041691| 0.0931691|
|ECDH10B_0788 | 1268.35915|       2.081515| 0.4589754| 3.260639| 0.0011116| 0.0320926|
|ECDH10B_0789 |  765.20941|       2.649082| 0.4304883| 4.794833| 0.0000016| 0.0000919|
|ECDH10B_0790 | 4079.52529|       1.885047| 0.3563550| 3.648286| 0.0002640| 0.0089146|

2. **genes-results.tsv**: a result table from a **DESeq** analysis giving *base means* across samples, *log2 fold changes*, *standard errors*, *test statistics*, *p-values* and *adjusted p-values*. The full table is available here [genes-results.tsv](https://github.com/noelnamai/dgexp/blob/master/results/genes-results.tsv).

|genes          | baseMean| log2FoldChange|    lfcSE|       stat|    pvalue| padj|
|:--------------|--------:|--------------:|--------:|----------:|---------:|----:|
|EBG00001089504 | 30.12564|     -1.0299979| 1.563530| -0.2846351| 0.7759237|    1|
|EBG00001089505 |  0.00000|             NA|       NA|         NA|        NA|   NA|
|EBG00001089506 |  0.00000|             NA|       NA|         NA|        NA|   NA|
|EBG00001089507 |  0.00000|             NA|       NA|         NA|        NA|   NA|
|EBG00001089508 | 22.26069|      1.1132234| 1.496692|  0.3529524| 0.7241241|    1|
|EBG00001089509 | 35.45083|     -0.3843498| 1.404188|  0.0000000| 1.0000000|    1|

3. **[ma-plot.png](https://github.com/noelnamai/dgexp/blob/master/results/ma-plot.png)**: the MA-plot from base means and log fold changes.

![Alt text](./results/ma-plot.png)

4. **[max-counts-plot.png](https://github.com/noelnamai/dgexp/blob/master/results/max-counts-plot.png)**: plot of normalized counts for a single gene with the maximum *padj*.

![Alt text](./results/max-counts-plot.png)

5. **[min-counts-plot.png](https://github.com/noelnamai/dgexp/blob/master/results/min-counts-plot.png)**: plot of normalized counts for a single gene with the minimum *padj*.

![Alt text](./results/min-counts-plot.png)

