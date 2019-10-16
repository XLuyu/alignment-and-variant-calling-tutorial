# NGS tutorial

# NGS alignment and variant calling

This tutorial is forked from [ekg](https://github.com/ekg/alignment-and-variant-calling-tutorial). I specialized it for computer science people who is new to NGS data analysis. This tutorial steps through some basic tasks in alignment and variant calling using a handful of Illumina sequencing data sets. I removed some trivial details for you to go through quickly.

## Part 0: Setup

We're going to use a bunch of fun tools for working with genomic data on Linux (For SoC student, you can use xcnc cluster). In fact, you don't really need to install all of them now at one go. Instead, You may want to install each of them when it is required.

1. [sra-tools](https://github.com/ncbi/sra-tools/wiki): Download publicly available dataset from NCBI
1. [bwa](https://sourceforge.net/projects/bio-bwa/files/): Align reads to reference genome
1. [samtools](https://sourceforge.net/projects/samtools/files/): Operate aligned read file.
1. [freebayes](https://github.com/ekg/freebayes/releases): Call variant
2. [IGV](http://software.broadinstitute.org/software/igv/download)

## Part 1: Aligning E. Coli data with `bwa mem`

[E. Coli K12](https://en.wikipedia.org/wiki/Escherichia_coli#Model_organism) is a common laboratory strain that has lost its ability to live in the human intestine, but is ideal for manipulation in a controlled setting.
The genome is relatively short, and so it's a good place to start learning about alignment and variant calling.

### E. Coli K12 reference genome

We'll get some test data to play with. First, [the E. Coli K12 reference genome](http://www.ncbi.nlm.nih.gov/nuccore/556503834), from NCBI. It's a bit of a pain to pull out of the web interface so [you can also download it here](http://hypervolu.me/~erik/genomes/E.coli_K12_MG1655.fa).

```bash
# the start of the genome, which is circular but must be represented linearly in FASTA
curl -s http://hypervolu.me/~Eerik/genomes/E.coli_K12_MG1655.fa | head
# >NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
# AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
# TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
# ...
```

As you see, genome is actually a string of A,C,T, and G (sometimes N for unknown). Currently, sequencing technologies are able to **randomly** recognize a segment on it (viz. a substring), with length of 100~300. People usually generate such substring million times so that genome is tiled 30 times in the dataset. Each substring is called **read**. 

For example, the picture below is a part of genome. You can see its sequence at the bottom. The grey bars in the middle are reads.

![](https://i.imgur.com/cfQyl6W.png)

### E. Coli K12 Illumina 2x300bp MiSeq sequencing results

For testing alignment, let's get some data from a [recently-submitted sequencing run on a K12 strain from the University of Exeter](http://www.ncbi.nlm.nih.gov/sra/?term=SRR1770413). To download it from NCBI, we need to use an executable `fastq-dump` in the `sra-tools` (Part 0, 1.). It allows directly downloading data from a particular sequencing run ID:

```bash
fastq-dump --split-files SRR1770413
```

The downloaded data files are in [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format) for further processing. The `--split-files` part of the command ensures we get two files, one for the first and second mate in each pair. We'll use them in this format when aligning.

### E. Coli O104:H4 HiSeq 2000 2x100bp

As a point of comparison, let's also pick up a [sequencing data set from a different E. Coli strain](http://www.ncbi.nlm.nih.gov/sra/SRX095630[accn]). This one is [famous for its role in foodborne illness](https://en.wikipedia.org/wiki/Escherichia_coli_O104%3AH4#Infection) and is of medical interest.

```bash
fastq-dump --split-files SRR341549
```

### Setting up our reference indexes

#### FASTA file index

First, we'll want to allow tools (such as our variant caller) to quickly access certain regions in the reference. This is done using the `samtools (Part 0, 3.)` .fai FASTA index format, which records the lengths of the various sequences in the reference and their offsets from the beginning of the file.

```bash
samtools faidx E.coli_K12_MG1655.fa
```

Now it's possible to quickly obtain any part of the E. Coli K12 reference sequence. For instance, we can get the 200bp from position 1000000 to 1000200. We'll use a special format to describe the target region `[chr]:[start]-[end]`.

```bash
samtools faidx E.coli_K12_MG1655.fa NC_000913.3:1000000-1000200
```

We get back a small FASTA-format file describing the region:

```text
>NC_000913.3:1000000-1000200
GTGTCAGCTTTCGTGGTGTGCAGCTGGCGTCAGATGACAACATGCTGCCAGACAGCCTGA
AAGGGTTTGCGCCTGTGGTGCGTGGTATCGCCAAAAGCAATGCCCAGATAACGATTAAGC
AAAATGGTTACACCATTTACCAAACTTATGTATCGCCTGGTGCTTTTGAAATTAGTGATC
TCTATTCCACGTCGTCGAGCG
```

#### BWA's FM-index

`BWA (Part 0, 2.)` uses the [FM-index](https://en.wikipedia.org/wiki/FM-index), which a compressed full-text substring index based around the [Burrows-Wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform).
To use this index, we first need to build it:

```bash
bwa index E.coli_K12_MG1655.fa
```

You should see `bwa` generate some information about the build process:

```text
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 2.26 seconds elapse.
[bwa_index] Update BWT... 0.04 sec
[bwa_index] Pack forward-only FASTA... 0.03 sec
[bwa_index] Construct SA from BWT and Occ... 0.72 sec
[main] Version: 0.7.8-r455
[main] CMD: bwa index E.coli_K12_MG1655.fa
[main] Real time: 3.204 sec; CPU: 3.121 sec
```

And, you should notice a new index file which has been made using the FASTA file name as prefix:

```bash
ls -rt1 E.coli_K12_MG1655.fa*
# -->
E.coli_K12_MG1655.fa
E.coli_K12_MG1655.fa.fai
E.coli_K12_MG1655.fa.bwt
E.coli_K12_MG1655.fa.pac
E.coli_K12_MG1655.fa.ann
E.coli_K12_MG1655.fa.amb
E.coli_K12_MG1655.fa.sa
```

### Aligning our data against the E. Coli K12 reference

Here's an outline of the steps we'll follow to align our K12 strain against the K12 reference:

1. use `bwa` to generate SAM records for each read
2. convert the output to BAM
3. sort the output

We could the steps one-by-one, generating an intermediate file for each step.
However, this isn't really necessary unless we want to debug the process, and it will make a lot of excess files which will do nothing but confuse us when we come to work with the data later.
Thankfully, it's easy to use [unix pipes](https://en.wikiepdia.org/wiki/Pipeline_%28Unix%29) to stream most of these tools together (see this [nice thread about piping bwa and samtools together on biostar](https://www.biostars.org/p/43677/) for a discussion of the benefits and possible drawbacks of this).

You can now run the alignment using a piped approach. **Replace `$threads` with the number of CPUs you would like to use for alignment.** Not all steps in `bwa` run in parallel, but the alignment, which is the most time-consuming step, does. You'll need to set this given the available resources you have.

```bash
bwa mem -t $threads \
    -R '@RG\tID:K12\tSM:K12' \
    E.coli_K12_MG1655.fa \
    SRR1770413_1.fastq.gz \
    SRR1770413_2.fastq.gz \
| samtools sort -@ $threads - -o SRR1770413.bam
```

Breaking it down by line:

- *alignment with bwa*: `bwa mem -t $threads -R '@RG\tID:K12\tSM:K12'` --- this says "align using so many threads" and also "label the reads as the read group K12 and the sample name K12"
- *reference and FASTQs* `E.coli_K12_MG1655.fa SRR1770413_1.fastq.gz SRR1770413_2.fastq.gz` --- this just specifies the base reference file name (`bwa` finds the indexes using this) and the input alignment files. The first file should contain the first mate, the second file the second mate.
- *sorting and conversion to BAM format*: `samtools sort -@ $threads -` --- this reads SAM from stdin (the `-` specifier in place of the file name indicates this) and sort reads using so many threads (by `-@` option). `-o SRR1770413.bam` specifies the output name and format (in BAM format).


Now, run the same alignment process for the O104:H4 strain's data. Make sure to specify a different sample name via the `-R '@RG...` flag incantation to specify the identity of the data in the BAM file header and in the alignment records themselves:

```bash
bwa mem -t $threads \
    -R '@RG\tID:O104_H4\tSM:O104_H4' \
    E.coli_K12_MG1655.fa \
    SRR341549_1.fastq.gz \
    SRR341549_2.fastq.gz \
| samtools sort -@ $threads - -o SRR341549.bam
```

As a standard post-processing step, it's helpful to add a BAM index to the files. This let's us jump around in them quickly using BAM compatible tools that can read the index. We use `samtools` to achieve this:

```bash
samtools index SRR1770413.bam
samtools index SRR341549.bam
```

## Part 2: Calling variants

Now that we have our alignments sorted, we can quickly determine variation against the reference by scanning through them using a variant caller.
There are many options, including [samtools mpileup](http://samtools.sourceforge.net/samtools.shtml), [platypus](http://www.well.ox.ac.uk/platypus), and the [GATK](https://www.broadinstitute.org/gatk/). For this tutorial, we'll keep things simple and use [freebayes](https://github.com/ekg/freebayes). 

It has a number of advantages in this context (bacterial genomes), such as long-term support for haploid (and polyploid) genomes. However, the best reason to use it is that it's very easy to set up and run, and it produces a very well-annotated VCF output that is suitable for immediate downstream filtering.

> [Comment] Note the original author of this tutorial is exactly the author of freebayes, ekg. As he said, freebayes is indeed the easiest one to start with. Only one-liner command is required. In Contrast, GATK requires a config file which may include hundreds of lines. But, GATK is the most commonly used one and recognized as the best.


### Variant calls with `freebayes`

It's quite easy to use `freebayes (Part 0, 4.)` provided you have your BAM file completed. We use `--ploidy 1` to indicate that the sample should be genotyped as haploid.

```bash
freebayes -f E.coli_K12_MG1655.fa --ploidy 1 SRR1770413.bam >SRR1770413.vcf
```
## Part 3: visualize your data

We use `IGV (Part 0, 5.)` to visualize all files in your hand, e.g. fasta, bam, and vcf. The installation is quite simple. But it requires X11 to show a window on your host, if you run it on a server without monitor. For windows user, it can be easily solved by using [`MobaXterm`](https://mobaxterm.mobatek.net/download.html) as terminal.

After IGV starts, you can see a window. Then,
1. load genome fasta file ![](https://i.imgur.com/ztm4SOB.png)
2. load bam and vcf file ![](https://i.imgur.com/9O9RZsz.png)

Play with it to get familiar. A SNP looks like this:

![](https://i.imgur.com/95Ys3rE.png)

While reference genome is T, this individual has 50% A and 50% T. This could make the guy has different hair color or something else.

## Part 4: Resource

* The human reference genome [hg38](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz). Note that E. coli has only one chromosome, but human has 22+XY. Each chromosome is an independent long string.

* To extract read information from BAM files in coding manner, an easy way is to learn Python and use the package [Pysam](https://pysam.readthedocs.io/en/latest/index.html)

