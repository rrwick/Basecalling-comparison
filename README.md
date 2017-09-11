<p align="center"><img src="images/logo.png" alt="logo" width="100%"></p>

This repo contains a comparison of available basecallers for Oxford Nanopore Technologies (ONT) sequencers. Basecallers, for those not familiar, are the programs which translate the raw electrical signal from an ONT sequencer to a DNA sequence. Basecalling is interesting because it's a hard problem (modern basecallers all seem to tackle it with neural networks) and it's a huge part of what makes ONT sequencing good or bad. Getting a piece of DNA through a pore and measuring the current is only half the battle; the other half is in the computer. It's a very active field, with both ONT themselves and other researchers developing methods.

For each basecaller, I assess the accuracy of the reads and the accuracy of an assembly. Read accuracy is interesting for obvious reasons – more accurate reads are nice! Assembly accuracy is also interesting because it provides a window into the nature of the basecalling errors. For example, consider a hypothetical set of reads with a mediocre accuracy of 85% but a truly random error profile (i.e. no systematic error). Despite their error rate, these reads could result in a perfect assembly because their errors are all 'averaged out' in the consensus. In contrast, now consider a set of reads with an excellent 98% accuracy but they all make the _same mistakes_ (i.e. error is all systematic, not random). The resulting assembly will also have a 98% error rate. Which read set is better? That probably depends on how you're using them, but in my line of work (assembly), I'd much prefer the first.

In particular, I hope these results help to answer the question: _Should I go back to old reads and re-basecall them?_ If basecallers improve enough, it might be worth it for the better reads. But that would take a lot of CPU time, so you probably don't want to do it unless there's a real benefit to be had.

The included scripts aren't intended for general use (i.e. this isn't a clone-it-and-run-it kind of repo), so if you want to reproduce my results you'll need to modify some paths in [`basecalling_comparison.sh`](basecalling_comparison.sh). I used an R9.4 dataset of _Klebsiella pneumoniae_ reads, so my results may be biased toward bacterial genomics of Enterobacteriaceae.

If you'd like to try this analysis using the same data, here are the relevant links:
* [Reference hybrid assembly](https://ndownloader.figshare.com/files/8810704)
* [Raw fast5 files](https://ndownloader.figshare.com/files/9199063)
* [Relevant repo/paper](https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing)



# Basecallers tested

For each basecaller I have only used the training models included with the program. Custom training of the neural net is out of scope for this analysis.


### Nanonet

[Nanonet](https://github.com/nanoporetech/nanonet) is ONT's first generation neural network basecaller. I used the most [latest commit](https://github.com/nanoporetech/nanonet/commit/a5a832bb3c82fbde091554142fab491d4bec2664) (at the time of writing) on the master branch. This seems to be functionally identical to v2.0.0, so that's what I've called it below. Nanonet no longer appears to be under active development, so this may be the last version.

```
nanonetcall --chemistry r9.4 --write_events --min_len 1 --max_len 1000000 --jobs 40 raw_fast5_dir > /dev/null
```
The `--min_len` and `--max_len` options were set so Nanonet wouldn't skip any reads. While Nanonet outputs its basecalls to stdout in fasta format, I've ignored that and instead used the `--write_events` options so it stores fastq basecalls in the fast5 files, which I can extract later. Unlike Albacore, which makes a copy of fast5 files, Nanonet modifies the original ones in the `raw_fast5_dir` directory.


### Albacore

Albacore is ONT's official command-line basecaller. I tested versions 0.8.4, 0.9.1, 1.0.4, 1.1.2, 1.2.6 and 2.0.1. The transducer basecaller (helps with homopolymers) was added in v1.0. Event-free basecalling first appears in v2.0. I think you need an account with the [Nanopore community](https://community.nanoporetech.com/) to get the download links for Albacore.

```
# Albacore v0.8.4 and v0.9.1:
read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i raw_fast5_dir -t 40 -s output_dir

# Albacore v1.0.4:
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5_dir -t 40 -s output_dir

# Albacore v1.1.2 and v1.2.6:
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5_dir -t 40 -s output_dir -o fast5

# Albacore v2.0.2:
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5_dir -t 40 -s output_dir -o fast5 --disable_filtering
```

Albacore v1.1 and later can basecall directly to fastq file with `-o fastq`, which saves disk space and is usually more convenient. However, for this experiment I used `-o fast5` to keep my analysis options open.


### Scrappie

[Scrappie](https://github.com/nanoporetech/scrappie) is ONT's research basecaller. Successful developments here seem to eventually work their way into Albacore. I tested versions 1.0.0 and 1.1.0.

Scrappie can be run as `scrappie events` (where it basecalls from event segmentation) or as `scrappie raw` (where it basecalls directly from the raw signal). For Scrappie v1.0.0, running as `scrappie events` relies on preexisting event data in the fast5s. For that test I used the fast5s as produced by Albacore 1.2.6 – the last Albacore version to do event segmentation.

In Scrappie v1.1.0, there are three different raw basecalling models to choose from (raw_r94, rgr_r94 and rgrgr_r94) and I tried each.

Unlike other basecallers, Scrappie does not produce fastq output, either directly or by writing it into the fast5 files.

```
# Scrappie v1.0.0:
scrappie events --albacore --threads 40 albacore_v1.2.6_fast5 > scrappie_v1.0.0_events.fasta
scrappie raw --threads 40 raw_fast5 > scrappie_v1.0.0_raw.fasta

# Scrappie v1.1.0:
scrappie events --threads 40 raw_fast5 > scrappie_v1.1.0_events.fasta
scrappie raw --model raw_r94 --threads 40 raw_fast5 > scrappie_v1.1.0_raw_raw_r94.fasta
scrappie raw --model rgr_r94 --threads 40 raw_fast5 > scrappie_v1.1.0_raw_rgr_r94.fasta
scrappie raw --model rgrgr_r94 --threads 40 raw_fast5 > scrappie_v1.1.0_raw_rgrgr_r94.fasta
```


### Chiron

[Chiron](https://github.com/haotianteng/chiron) is a third-party neural-network basecaller. The first release (v0.1) did not work on my reads (`extract_sig_ref.py` failed to get the signal from the fast5 files) so I instead used the [latest commit](https://github.com/haotianteng/chiron/commit/847ad101d338a253a7dfed2023c3ed8758886e7e) (at the time of writing) on the master branch.

```
chiron call -i raw_fast5 -o output_dir
```


### Excluded basecallers

I did not include [basecRAWller](https://basecrawller.lbl.gov/) for two reasons. First, as stated in [its paper](http://www.biorxiv.org/content/early/2017/05/01/133058), basecRAWller's has a somewhat different focus from the other programs: streaming basecalling. Second, it is not freely available to downloaded without registration on the website.

I _am_ interested in [Guppy](https://github.com/nanoporetech/guppy), an ONT basecaller designed for fast GPU-accelerated performance. However, it is only available to users who have signed the ONT developer license agreement (which is why you might have gotten a 404 if you just tried that link). If and when Guppy becomes publicly available, I'll add it to this comparison.



# Method

### Basecalling

My sample data is from a barcoded run, so I first collected the fast5 files which I had previously binned to the appropriate sample. I also tossed out any fast5 files less than 100 kb in size – this was to subsample for longer reads in a manner that's hopefully not biased towards any particular basecaller.

The commands to do basecalling vary depending on the program. Where possible I basecalled to fast5 files so I could later use Nanopolish. Whenever low-level tuning parameters were available, I stuck with the defaults. To extract a fastq from the fast5s, I used my own script [`fast5_to_fastq.py`](https://github.com/rrwick/Fast5-to-Fastq), but many other tools exist to do the same job.


### Read identity

To assess read identity, I aligned the reads to the reference using [minimap2](https://github.com/lh3/minimap2) and then used [`read_length_identity.py`](read_length_identity.py) to get a per-read summary. Only the aligned parts of the read were used to calculate the read's identity. If less than 50% of the read aligned, it was deemed unaligned. This script also determines the read length to reference length ratio for each read, to see if insertions or deletions are more likely. 


### Assembly identity

Before assembly, I used [Porechop](https://github.com/rrwick/Porechop) and [Filtlong](https://github.com/rrwick/Filtlong) to clean up the read set*. I then assembled with [Unicycler](https://github.com/rrwick/Unicycler). Unicycler conducts multiple rounds of [Racon](https://github.com/isovic/racon), so the final assembly accuracy is defined by the Racon consensus (which in my experience is a bit higher accuracy than a [Canu](https://github.com/marbl/canu) assembly).

To get a distribution of assembly identity, I used [`chop_up_assembly.py`](chop_up_assembly.py) to divide the assembly into 10 kbp 'reads' which were assessed in the same way as the actual reads ([`read_length_identity.py`](read_length_identity.py)). When possible, I then used [Nanopolish](https://github.com/jts/nanopolish) [v0.7.1](https://github.com/jts/nanopolish/releases/tag/v0.7.1) to improve the assembly and assessed identity again.

<sup>* I did use Filtlong's reference-based mode (using Illumina reads) with trimming and splitting, so this assembly method isn't _truly_ ONT-only. However, I found that doing so led to more consistent assemblies (always getting one contig per replicon) which made it a lot easier to compare them.</sup>


# Results


### Total yield

_YIELD BAR PLOT HERE_

You might expect that each basecaller would produce approximately the same total yield. E.g. a read that's 10 kbp in one basecaller would be about 10 kbp in each basecaller. That's true for the most part, but Nanonet is a notable exception. For most reads, it produced a much shorter sequence than other basecallers, sometimes drastically so. For example, all versions of Albacore basecalled one read (`d2e65643-98b4-4a61-ad36-119e55250b28`) to a 34+ kbp sequence. Nanonet produced a 518 bp sequence for the same read. I don't have an explanation for this odd behaviour.

Other oddities you might notice are Albacore v0.9.1, which produced noticably more sequence than the other basecallers, and Scrappie events, which produced less. That will explained in the 'Relative read length' section.


### Read identity

_VIOLIN PLOT HERE_

This first analysis tackles the most obvious question: how accurate are the basecalled reads? The plot above shows the read identity distribution, with the median (weighted by read length) marked as a horizontal line. Reads that could not be aligned were given an identity of 0%, creating a bulge at the bottom of the distribution. Reads with an actual identity below about 65% usually to fail to align and fall to 0.

Nanonet performed poorly, with a low median identity and a significant proportion of unaligned reads. Its curiously high peak of about 99% is misleading and results from its short output sequences discussed above. While a few Nanonet 'reads' did indeed align to the reference with up to 99% identity, these were really small fragments (hundreds of bp) of larger reads.

Albacore v0.9.1 is quite peculiar - its read identity distribution has a wide variance and the lowest median of all tested basecallers. All other versions of Albacore performed well, and reassuringly the latest version (v2.0.2) performed best. But surprisingly, v0.8.4 performed second best! The transducer algorithm introduced in v1 did not seem to create a large benefit, possibly because my data isn't too homopolymer-heavy.

Scrappie also did fairly well, overall comparible in accuracy to Albacore.
_WRITE MORE HERE WHEN SCRAPPIE RAW v1.0.0 IS DONE_
Interestingly, Scrappie produced significantly more unalign-able reads than Albacore (excluding v0.9.1). This is true even for Scrappie raw v1.1.0 rgr_r94, which had the best overall median identity. It 



### Relative read length

_VIOLIN PLOT HERE_

By looking at the relative length of read to reference in each alignment, we can get a picture of whether the basecaller is more prone to insertions or deletions. 100% (same length) means that insertions and deletions are equally likely. <100% means that deletions are more common than insertions. >100% means that insertions are more common than deletions. Curiously, many basecallers produce a distinctly bimodal distribution here. For example, Albacore v2.0.2 has a good median value of quite near 100%, but any individual read is more likely to fall either around 99% or 101%. I'm not sure what might be behind this.

Once again, Albacore v0.9.1 stands out as unusual, tending to overly-long reads, but in a bimodal way. Scrappie events tends to short reads. These explain the total yield differences we saw earlier.



### Assembly identity

_IDENTITY VIOLIN PLOT HERE_

This analysis is my personal favourite: how accurate are the _consensus_ sequences? I don't particularly care if individual reads have low identity if they produce an accurate assembly.

The first thing you might notice is that the difference between basecallers is more pronounced here than it was for read identity. Most surprisingly, Albacore v0.9.1, which was the worst performing in read identity, produces very good assemblies! This suggests that while its reads have a high degree of error, they have relatively little _systematic error_, and so the errors tend to average out in the consensus. Another notable finding is that while Scrappie events and Scrappie raw are comparible in read identity, Scrappie raw is clearly superior for assembly identity.

_WRITE MORE HERE WHEN SCRAPPIE RAW v1.0.0 IS DONE_

To better see how read identity and assembly identity relate, here are the median identities in a scatter plot:

__SCATTER PLOT HERE__

If the two values were tightly correlated, the points would fall on a diagonal line. They mostly do, with Albacore v0.9.1 being the most obvious exception.

Finally, it's interesting to look at the assembly insertions-vs-deletions biases, like we did for reads:

_RELATIVE LENGTH VIOLIN PLOT HERE_



### Post-Nanopolish assembly identity

_VIOLIN PLOT HERE_



# Conclusions

### Recommendations

For most users, my current recommendation would be to just use the latest version of Albacore: v2.0.2. It does well with both read and assembly accuracy. Scrappie raw is also a good choice, especially with the rgr_r94 model, though it doesn't offer much advantage over Albacore.

The recommendation would have been harder before Albacore v2.0.2 was released. Then, the best assembly accuracy would be produced by the basecaller with the _worst_ read accuracy: Albacore v0.9.1. This apparent paradox can be explained by that version having lots of random error but relatively low systematic error. Whether or not it would be good choice might depend on the intricacies of your analysis. My hunch is that it may often come down to read depth: low-depth analyses benefit from read accuracy while high-depth analyses benefit from consensus accuracy. We've dodged that tough decision for the moment (just use v2.0.2), but may be faced with a similar dilemma if a future basecaller excels at consensus accuracy over read accuracy or vice versa.

Finally, Nanonet seems a bit dated and should probably be avoided. However, it does allow for custom training of the neural network and may therefore be of interest to power users interested in experimenting.


### Future work

_My_ future work is easy: trying new basecallers as they are released and adding them to this repo. In particular, I'm interested to try Guppy, ONT's GPU basecaller, when that is made public. Check back occasionally for new data!

The much harder task lies with the basecaller authors: reducing systematic error. As it currently stands, systematic basecalling errors lead to residual errors in assemblies, even after Nanopolish. This makes it hard to recommend an ONT-only approach for many types of genomics where accuracy matters (read more in our paper on this topic). If systematic error can be eliminated, then ONT-only assemblies will approach 100% accuracy, and then ONT will be true Illumina alternative.

