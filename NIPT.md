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

## nucleosome
DNA wraps around nucleosome, which prevents DNA from degradation. Fetal DNA suffers more from degradation and thus reads from fetus start a little bit closer to center of nucleosome. By inferring positions of nucleosome, a distribution of read start point (distance to the center of the nearest nucleosome) can be obtained and used as feature.

## AI potential
Nowadays, AI usually refers to deep learning. It provides better modeling ability. I have tried a naive DNN to perform FF prediction and yielded a slightly better result than Linear/Ridge regression. However, the main challenge in this project is feature engineering. We only have hundreds of samples, but read start point distribution also ranges from -150 to 150. This means feature is in the same magnitude of samples. But machine learning is practical when sample size is much larger then feature size.

## Have a glance at data
To get familar with NGS data (Actually quite simple), you may want to go through a short tutorial https://github.com/taoistly/alignment-and-variant-calling-tutorial . For NIPT, part 2 and later sections are not very relevant.
