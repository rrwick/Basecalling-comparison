<p align="center"><img src="images/logo.png" alt="logo" width="100%"></p>

This experiment/repository/blog post aims to compare available basecallers for Oxford Nanopore data. In particular, I'm interested in the identity of the reads and of the resulting assembly. I used an R9.4 dataset for _Klebsiella pneumoniae_, so my results may be biased toward bacterial genomics of Enterobacteriaceae.


# Basecallers tested

### Albacore

I tried a range of versions: 0.8.4, 0.9.1, 1.0.4, 1.1.2, 1.2.6, 2.0.0

The transducer basecaller which helps with homopolymers was added in v1.0. Event-free basecalling first appears in v2.0. 

I'm not aware of publicly available download links for Albacore - I think you need an account with the [Nanopore community](https://community.nanoporetech.com/) to get it.


### Scrappie

[Scrappie](https://github.com/nanoporetech/scrappie) is ONT's research basecaller. Successful developments here seem to eventually work their way into Albacore.

I used the [latest commit](https://github.com/nanoporetech/scrappie/commit/16b9f366694689cd51ba1c134444bccc31aa4b80) on the master branch at the time of writing this.


### Nanonet

[Nanonet](https://github.com/nanoporetech/nanonet) is ONT's first generation neural network basecaller.

I used the [latest commit](https://github.com/nanoporetech/nanonet/commit/a5a832bb3c82fbde091554142fab491d4bec2664) on the master branch at the time of writing this.


### Chiron

[Chiron](https://github.com/haotianteng/chiron) is a new third-party neural-network bascaller. I used the [latest commit](https://github.com/haotianteng/chiron/commit/5914d6cdfa59e8b7299c55ca9ed4d0115d3a3626) on the master branch at the time of writing this.



# Method

The script [`basecalling_comparison.sh`](basecalling_comparison.sh) contains the steps used in this analysis. Near the top of the script are some paths specific to my server, so if you want to reproduce the analysis on your own data, you'll need to modify those. I used the barcode 1 reads from [this data set](https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing).


### Collect fast5 files

My sample data is from a barcoded run, so I first collected just the fast5 files which were confidently binned to the appropriate sample. I then tossed out any fast5 files less than 100 kb in size - this was for the sake of subsampling the reads in a manner that's hopefully not biased towards any particular basecaller.


### Basecalling

The commands to do basecalling vary depending on the program. Where possible I basecalled to fast5 files so I can later use Nanopolish.

To extract a fastq from the fast5s, I used my own script [`fast5_to_fastq.py`](https://github.com/rrwick/Fast5-to-Fastq), but many tools exist to do the same job.


### Read identity

To assess read identity, I aligned the reads using [minimap2](https://github.com/lh3/minimap2) and then used [`read_length_identity.py`](read_length_identity.py) to get a per-read summary. If less than 50% of the read aligned, that script deems the read unaligned. If 50% or more aligned, then it uses only the aligned parts to get an average identity. This script also determines the read length to reference length ratio per read, to see if insertions or deletions are more likely. 


### Assembly identity

Before assembly, I used [Porechop](https://github.com/rrwick/Porechop) and [Filtlong](https://github.com/rrwick/Filtlong) to clean up the read set. I then assembled with [Unicycler](https://github.com/rrwick/Unicycler). Unicycler conducts multiple rounds of [Racon](https://github.com/isovic/racon), so the final assembly accuracy is defined by the Racon consensus (which in my experience is a bit better than an unpolished [Canu](https://github.com/marbl/canu) assembly).

To get a distribution of assembly identity, I used [`chop_up_assembly.py`](chop_up_assembly.py) to divide the assembly into 10 kbp 'reads' which were assessed in the same way as the actual reads ([`read_length_identity.py`](read_length_identity.py)).

When possible, I then used [Nanopolish](https://github.com/jts/nanopolish) [v0.7.1](https://github.com/jts/nanopolish/releases/tag/v0.7.1) to improve the assembly and assessed identity again. This was only possible with some of the basecallers.


# Results

### Read identity


### Read length


### Assembly identity


### Post-Nanopolish assembly identity





