##R.A.G.E. : Ribosome Assisted Genome Exploration

###Requirements

* Linux Operating System (64-bit)
* Python version 2.7 or higher
* BioPython module
* Cutadapt-1.4.2 or higher (for removing the Adaptors and cleaning the reads)
* Bowtie-1.0.1 or higher - (for aligning the reads to the genome)
* Samtools-0.1.19 or higher - (for converting sam files to bam files)
* CEM -v0.9.1 or higher - (for assembling transcripts from RNA-Seq reads)
* glimmer3.02 - (for de novo gene prediction)
* cufflinks-2.2.0 - (for converting gtf to gff3 file)

###Introduction

With the development of high-throughput sequencing technologies in past years, ribosome profiling has emerged as a novel technique that provides a quantitative and qualitative measure of mRNA translation and protein production. Ribosome profiling (or Ribo-Seq) is performed via nuclease footprinting, which involves digestion of all RNA segments unoccupied by the ribosomes. The ribosome protected fragments or “footprints” are purified and then subjected to deep sequencing. Given that this technique is rather new, few genome annotation pipelines have incorporated the information that Ribo-Seq provides into their predictions. RAGE is an automated software pipeline targeted towards prokyaryotic genomes that uses mRNA and Ribo-seq data along with a reference genome to produce a set of gene predictions that are ranked according to which data set aided the gene discovery. This pipeline was original created by the contributors for the course 03-713 at Carnegie Mellon University.

###Setup

If you haven’t already setup the paths for the tools, run the command below to input the paths directly into the bash_profile file.

###Example Commands
```
nohup python rage-v2.4.py $HOME/data/ecolik12/k12rnaseq.fastq $HOME/data/ecolik12/k12riboseq.fastq $HOME/data/ecolik12/k12genome.fasta &
```
```
nohup python rage-v2.4.py $HOME/data/ecolik12/k12rnaseq.fastq $HOME/data/ecolik12/k12riboseq.fastq $HOME/data/ecolik12/k12genome.fasta -d $HOME/results/data/outputsfinaltest/ -w 9#3#1 -s ecolik12 &
```

```
nohup python rage-v2.4.py $HOMEresults/2014-04-07/rnaseqtrimmed.fastq $HOME/results/2014-04-16/outputs/cutadapt/k12riboseq-ribotrimmed.fastq $HOME/data/ecolik12/k12genome.fasta -c false &
```
```
nohup python rage-v2.4.py $HOMEresults/2014-04-07/rnaseqtrimmed.fastq $HOME/results/data/outputsfull/filter/ribotrimmedfiltered_rRNA_tRNA.fastq $HOME/data/ecolik12/k12genome.fasta -c false -f false &
```
```
nohup python rage-v2.4.py $HOME/results/data/outputfullwglimmer/bowtie/rnaseqtrimmed.sorted.bam $HOME/results/data/outputsfull/rnaseqtrimmed.sam $HOME/data/ecolik12/k12genome.fasta -c false -f false -b false &
```
###Parameters
```
Usage: rage-v2.4.py [-h] [-c {true,false}] [-f {true,false}] [-b {true,false}]
                    [-d WORKDIRCTORY] [-a RNA_ADAPTER] [-i RIBOSOME_ADAPTER]
                    [-w WEIGHTS] [-s SPECIES]
                 
                    Positional Arguments:
  rnainput              input RNA-seq fastq file
  riboinput             input Ribosome fastq file
  genomeinput           input Genome fasta file

Optional Arguments:
  -h, --help            show this help message and exit
  -c {true,false}, --cutadapt {true,false}
                        whether use cutadapt, if choose false, please input
                        the clean RNA-seq and ribosome fastq file; default is
                        true
  -f {true,false}, --filter_rRNA_tRNA {true,false}
                        whether filter rRNA and tRNA, if choose false, please
			input the clean ribosome fastq file; default is true
  -b {true,false}, --bowtie {true,false}
                        whether use bowtie to map, if choose false, please
                        input sam file for Ribo-footprinting and sorted .bam
                        for RNA-seq; default is true
  -d WORKDIRCTORY, --WorkDirctory WORKDIRCTORY
                    rnainput riboinput genomeinput
  -a RNA_ADAPTER, --rna_adapter RNA_ADAPTER
                        input custom adapt sequence for RNA-seq
  -i RIBOSOME_ADAPTER, --ribosome_adapter RIBOSOME_ADAPTER
                        input custom adapt sequence for Ribosome profiling
  -w WEIGHTS, --weights WEIGHTS
                        Assign points for ribo-footprinting, RNAtranscripts
                        and naiveORFs; format XX#XX#XX
  -s SPECIES, --species SPECIES
                        Species name of the data
```

###References

####LICENSE
This program is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

####CONTRIBUTORS
* Dr. C. Joel McManus
* Minh Le
* Easwaran Ramamurthy
* Xinyu Miao
* Wenzhu Liu
* Siddharth Gurdasani
* Patrick Ropp
