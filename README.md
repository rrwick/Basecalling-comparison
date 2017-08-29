<p align="center"><img src="images/logo.png" alt="logo" width="100%"></p>

This repo is something of a blog plot and aims to compare available basecallers for Oxford Nanopore data. Basecallers, for those not familiar, are the programs which translate the raw electrical signal from an Oxford Nanopore sequencer to a sequence. It's not an easy task, and modern basecallers all seem to tackle this problem using neural networks.

For each basecaller, I assess the accuracy of the reads and the accuracy of an assembly. Read accuracy is interesting for obvious reasons – more accurate reads are nice! Assembly accuracy is also interesting because it provides a window into the nature of the basecalling errors. For example, consider a hypothetical set of reads with a mediocre accuracy of 85% but a truly random error profile (i.e. no systematic error). Despite their error rate, these reads could result in a perfect assembly because their errors are all 'averaged out' in the consensus. In contrast, now consider a set of reads with an excellent 98% accuracy but they all make the _same mistakes_ (i.e. error is all systematic, not random). The resulting assembly will also have a 98% error rate. Which read set is better? That probably depends on how you're using them, but in my line of work (assembly), I'd much prefer the first.

The included scripts aren't intended for general use (i.e. this isn't a clone-it-and-run-it kind of repo), so if you want to reproduce my results you'll need to modify some paths in [`basecalling_comparison.sh`](basecalling_comparison.sh). I used an R9.4 dataset of _Klebsiella pneumoniae_ reads, so my results may be biased toward bacterial genomics of Enterobacteriaceae.

If you'd like to try this analysis using the same data, here are the relevant links:
* [Reference hybrid assembly](https://ndownloader.figshare.com/files/8810704)
* [Raw fast5 files](https://ndownloader.figshare.com/files/)
* [Relevant repo/paper](https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing)



# Basecallers tested

### Albacore

I tried a range of versions: 0.8.4, 0.9.1, 1.0.4, 1.1.2, 1.2.6, 2.0.0, 2.0.1

The transducer basecaller (helps with homopolymers) was added in v1.0. Event-free basecalling first appears in v2.0. 

I'm not aware of publicly available download links for Albacore – I think you need an account with the [Nanopore community](https://community.nanoporetech.com/) to get them.


### Scrappie

[Scrappie](https://github.com/nanoporetech/scrappie) is ONT's research basecaller. Successful developments here seem to eventually work their way into Albacore.

I used the [latest commit](https://github.com/nanoporetech/scrappie/commit/16b9f366694689cd51ba1c134444bccc31aa4b80) on the master branch at the time of writing this.


### Nanonet

[Nanonet](https://github.com/nanoporetech/nanonet) is ONT's first generation neural network basecaller.

I used the [latest commit](https://github.com/nanoporetech/nanonet/commit/a5a832bb3c82fbde091554142fab491d4bec2664) on the master branch at the time of writing this.


### Chiron

[Chiron](https://github.com/haotianteng/chiron) is a new third-party neural-network basecaller. I used the [latest commit](https://github.com/haotianteng/chiron/commit/5914d6cdfa59e8b7299c55ca9ed4d0115d3a3626) on the master branch at the time of writing this.



# Method

### Basecalling

My sample data is from a barcoded run, so I first collected just the fast5 files which I had previously binned to the appropriate sample. I also tossed out any fast5 files less than 100 kb in size – this was to subsample for longer reads in a manner that's hopefully not biased towards any particular basecaller.

The commands to do basecalling vary depending on the program. Where possible I basecalled to fast5 files so I could later use Nanopolish. To extract a fastq from the fast5s, I used my own script [`fast5_to_fastq.py`](https://github.com/rrwick/Fast5-to-Fastq), but many other tools exist to do the same job.


### Read identity

To assess read identity, I aligned the reads to the reference using [minimap2](https://github.com/lh3/minimap2) and then used [`read_length_identity.py`](read_length_identity.py) to get a per-read summary. Only the aligned parts of the read were used to calculate the read's identity. If less than 50% of the read aligned, it was deemed unaligned. This script also determines the read length to reference length ratio for each read, to see if insertions or deletions are more likely. 


### Assembly identity

Before assembly, I used [Porechop](https://github.com/rrwick/Porechop) and [Filtlong](https://github.com/rrwick/Filtlong) to clean up the read set. I then assembled with [Unicycler](https://github.com/rrwick/Unicycler). Unicycler conducts multiple rounds of [Racon](https://github.com/isovic/racon), so the final assembly accuracy is defined by the Racon consensus (which in my experience is a bit higher accuracy than a [Canu](https://github.com/marbl/canu) assembly).

To get a distribution of assembly identity, I used [`chop_up_assembly.py`](chop_up_assembly.py) to divide the assembly into 10 kbp 'reads' which were assessed in the same way as the actual reads ([`read_length_identity.py`](read_length_identity.py)).

When possible, I then used [Nanopolish](https://github.com/jts/nanopolish) [v0.7.1](https://github.com/jts/nanopolish/releases/tag/v0.7.1) to improve the assembly and assessed identity again. This was only possible with some of the basecallers.



# Results

### Read identity


### Unaligned reads


### Relative read length

The absolute read lengths from one basecaller to another are quite similar, as you might expect. What I'm looking at here is the length of the read relative to its alignment on the reference. 100% (same length) means that insertions and deletions are equally likely. <100% means that deletions are more common than insertions. >100% means that insertions are more common than deletions.


### Assembly identity


### Post-Nanopolish assembly identity



# Conclusions

