# NIPT Project
## Background
* Human genome consists of 24 chromosome. Each is a string(50~200MB) on {A,C,T,G}.
* NGS (Next-Generation Sequencing) can read millions of small substrings (random position of chromosome, ~100 in length) from sample's genome.
* By mapping those reads to reference genome, we can know the origin of each read.
* When enough reads are collected, many applications are avaiable. NIPT is one of them.

## NIPT 
> Fact: In pregnant woman's boold, 5%~25% cfDNA come from fetus.

1. collect blood from pregnant woman
2. use NGS to get reads
3. map read to reference genome
4. use statistics to check if any chromosome gets uncommon read count:
    * too much (z-score > 3.0): trisomy detected
    * too less (z-score < -3.0): monosomy detected

`NIPT is not too hard, simple statistics can offer acceptable/reliable result.`

## Fetal Fraction
FF (Fetal Fraction) is the proportion of fetus DNA in mother's blood. Usually 5%~25%, depends on a lot of factors.

FF plays a role in:
* Quality Control: if a sample has FF < 5%, then a resample is required
* Trisomy Detection: Higher FF suggests higher threshold for trisomy read count

Compared to NIPT, FF detection is much harder. It has to exploit some feature which is different between mother and fetus.

Currently, I'm working on this, using NGS data to predict FF. 

The feature I exploit is the fact that: Due to DNA degradation, for most positions on chromosomes, mother and fetus has different probability to get a sequencing read. 

With hundreds of samples and millions of reads for each sample, machine learning may be helpful to find out which regions/positions show significant difference on sequencing probability. Then, those positions can be used to make FF prediction.

## Have a glance at data
To get familar with NGS data (Actually quite simple), you may want to go through a short tutorial https://github.com/taoistly/alignment-and-variant-calling-tutorial . For NIPT, part 2 and later sections are not very relevant.
