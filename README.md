# Comparison of Oxford Nanopore basecalling tools

__Ryan R. Wick, Louise M. Judd and Kathryn E. Holt__
<br>
<sub>Department of Biochemistry and Molecular Biology, Bio21 Molecular Science and Biotechnology Institute,
<br>
University of Melbourne, Parkville, Victoria 3010, Australia</sub>


<p align="center"><img src="images/logo.png" alt="logo" width="100%"></p>



This repo contains a comparison of available basecallers for Oxford Nanopore Technologies (ONT) sequencing reads. It's a bit like a mini-paper, but I decided to put it here on GitHub (instead of somewhere more papery like [bioRxiv](http://www.biorxiv.org/)) so I can come back and update it as new versions of basecallers are released.

Basecallers, for those not familiar, are the programs which translate the raw electrical signal from an ONT sequencer to a DNA sequence. Basecalling is interesting because it's a hard machine learning problem (modern basecallers all seem to tackle it with neural networks) and because it's a huge part of what makes ONT sequencing good or bad. Getting a piece of DNA through a pore and measuring the current is only half the battle; the other half is in the computer. It's a very active field, with both ONT themselves and independent researchers developing methods.

For each basecaller, I assess the accuracy of the reads and of the resulting assembly. Read accuracy is interesting for obvious reasons – more accurate reads are nice! Assembly accuracy is interesting because shows whether the read errors can 'average out' with depth. In doing so it provides a window into the nature of the basecalling errors. For example, consider a hypothetical set of reads with a mediocre accuracy of 85% but a truly random error profile (i.e. no systematic error). Despite their error rate, these reads could result in a perfect assembly. Now consider a set of reads with an excellent 98% accuracy but they all make the _same mistakes_ (i.e. error is all systematic, not random) – their assembly will also have a 98% error rate. Which read set is better? That probably depends on how you're using them, but in my line of work, I'd prefer the first.

In particular, I hope these results help to answer the question: _Should I go back to old reads and re-basecall them with a newer basecaller?_ Doing so would take a lot of CPU time, so you probably don't want to do it unless there's a significant improvement.

As a final note, I used an R9.4 1D dataset of _Klebsiella pneumoniae_ reads for this analysis, so my results may well be biased toward that kind of data. I'm not sure how consistent these results are with other data types, e.g. eukaryote genomes, R9.5 flow cells, 1D<sup>2</sup> kits, etc.






## Data availability

I've included the scripts I used for basecalling and analysis in this repo, but you'll need to modify some paths in [`basecalling_comparison.sh`](basecalling_comparison.sh) to use it for yourself.

If you'd like to try this analysis using the same data, here are the relevant links:
* [Reference hybrid assembly](https://ndownloader.figshare.com/files/8810704)
* [Raw fast5 files](https://ndownloader.figshare.com/files/9199063)
* [Relevant repo/paper](https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing)






## Basecallers tested

For each basecaller I have only used the training model(s) included with the program. Custom training of the neural net is out of scope for this analysis. Similarly, whenever low-level tuning parameters were available, I stuck with the defaults.



### Nanonet

[Nanonet](https://github.com/nanoporetech/nanonet) is ONT's first generation neural network basecaller. I used the [latest commit](https://github.com/nanoporetech/nanonet/commit/a5a832bb3c82fbde091554142fab491d4bec2664) (at the time of writing) on the master branch. This seems to be functionally identical to v2.0.0, so that's what I've called it. Nanonet no longer appears to be under active development, so this may be the final version.

```
nanonetcall --chemistry r9.4 --write_events --min_len 1 --max_len 1000000 --jobs 40 raw_fast5_dir > /dev/null
```
I set the `--min_len` and `--max_len` options so Nanonet wouldn't skip any reads. While Nanonet outputs its basecalls to stdout in fasta format, I've ignored that and instead used the `--write_events` options so it stores fastq basecalls in the fast5 files, which I can extract later. Unlike Albacore, which makes a copy of fast5 files, Nanonet modifies the original ones in the `raw_fast5_dir` directory.



### Albacore

Albacore is ONT's official command-line basecaller. I tested versions 0.8.4, 0.9.1, 1.0.4, 1.1.2, 1.2.6 and 2.0.2. All of these versions were released in 2017, which shows just how rapidly basecaller development is progressing. The transducer basecaller (helps with homopolymers) was added in v1.0. Basecalling from raw signal (without segmenting the signal into events) first appears in v2.0. Albacore can be downloaded from the [Nanopore community](https://community.nanoporetech.com/downloads), but you'll need an account to log in.

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

Albacore v1.1 and later can basecall directly to fastq file with `-o fastq`. This saves disk space and is usually more convenient (especially since Nanopolish v0.8), but for this experiment I used `-o fast5` to keep my analysis options open.



### Scrappie

[Scrappie](https://github.com/nanoporetech/scrappie) is ONT's research basecaller. Successful developments here seem to eventually work their way into Albacore. I tested versions 1.0.0 and 1.1.0.

Scrappie can be run as `scrappie events` (where it basecalls from event segmentation) or as `scrappie raw` (where it basecalls directly from the raw signal). For Scrappie v1.0.0, running as `scrappie events` relies on pre-existing event data in the fast5s, so I used the fast5s produced by Albacore 1.2.6 – the last Albacore version to do event segmentation. In Scrappie v1.1.0, there are three different raw basecalling models to choose from (raw_r94, rgr_r94 and rgrgr_r94) and I tried each.

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

Unlike other basecallers, Scrappie does not have fastq output, either directly or by writing it into the fast5 files. It only produces fasta reads.



### Excluded basecallers

I am currently trying two other basecallers: [basecRAWller](https://basecrawller.lbl.gov/) and [Chiron](https://github.com/haotianteng/chiron). I'll add their results after they finish. I am also interested in [Guppy](https://github.com/nanoporetech/guppy), an ONT basecaller designed for fast GPU-accelerated performance. However, it is only available to users who have signed the ONT developer license agreement (which is why you might have gotten a 404 if you just tried that link). If and when Guppy becomes publicly available, I'll add it to this comparison too.

Unfortunately, I cannot compare with the old cloud-based Metrichor basecalling, as it's no longer available. I also cannot test the integrated basecalling in MinKNOW (ONT's sequencing software). I believe MinKNOW's integrated basecalling shares much in common with Albacore, but I don't know which Albacore versions correspond to which MinKNOW versions.






## Method

### Sequencing

These reads are the same ones used in our recent paper: [Completing bacterial genome assemblies with multiplex MinION
sequencing](http://www.biorxiv.org/content/early/2017/07/07/160614). Look there if you're interested in the wet-lab side of things.



### Read preparation

I first collected the fast5 files which I had previously binned to the barcode 1. I also tossed out any fast5 files less than 100 kilobytes in size – this was to subsample for longer reads in a manner that's hopefully not biased towards any particular basecaller. To extract a fastq from the fast5s, I used my own script [`fast5_to_fastq.py`](https://github.com/rrwick/Fast5-to-Fastq), but many other tools exist to do the same job.



### Read identity

To assess read identity, I aligned the reads to the reference using [minimap2](https://github.com/lh3/minimap2) (lovely tool, very fast) and then used [`read_length_identity.py`](read_length_identity.py) to get a per-read summary. Only the aligned parts of the read were used to calculate the read's identity. The definition used for 'identity' is the same as how BLAST defines it: the number of matching bases in the alignment divided by the total bases in the alignment (including gaps). If less than 50% of the read aligned, it was deemed unaligned and given an identity of 0%. This script also determines the read length to reference length ratio for each read, to see if insertions or deletions are more likely.



### Assembly identity

Before assembly, I used [Porechop](https://github.com/rrwick/Porechop) and [Filtlong](https://github.com/rrwick/Filtlong) to clean up the read set\*. I then assembled with [Unicycler](https://github.com/rrwick/Unicycler), which conducts multiple rounds of [Racon](https://github.com/isovic/racon), so the final assembly accuracy is defined by the Racon consensus (which in my experience is a bit higher accuracy than a [Canu](https://github.com/marbl/canu) assembly). To get a distribution of assembly identity, I used [`chop_up_assembly.py`](chop_up_assembly.py) to divide the assembly into 10 kbp 'reads' which were assessed in the same way as the actual reads ([`read_length_identity.py`](read_length_identity.py)).

<sup>\* I did use Filtlong's reference-based mode (using Illumina reads) with trimming and splitting, so this assembly method isn't _truly_ ONT-only. However, I found that doing so led to more consistent assemblies (always getting one contig per replicon) which made it a lot easier to compare them.</sup>



### Nanopolish

I used [Nanopolish](https://github.com/jts/nanopolish) [v0.8.1](https://github.com/jts/nanopolish/releases/tag/v0.7.1) to improve each assembly and then assessed identity again. Since v0.8.0, Nanopolish can be run without event-data-containing fast5 files, which lets me run it with any basecaller! However, for non-Albacore basecallers I did have to alter read names so they were Albacore-like and compatible with the `nanopolish index` command.






## Results

### Total yield

<p align="center"><img src="images/total_yield.png" width="90%"></p>

You might expect that each basecaller would produce approximately the same total yield. E.g. a read that makes a 10 kbp sequence in one basecaller would be about 10 kbp in each basecaller. That's true for the most part, but Nanonet is a notable exception. For most reads, it produced a much shorter sequence than other basecallers, sometimes drastically so. For example, all versions of Albacore basecalled one read (d2e65643-98b4-4a61-ad36-119e55250b28) to a 34+ kbp sequence. Nanonet produced 518 bp for the same read. I don't have an explanation for this odd behaviour.

Other oddities you might notice are Albacore v0.9.1, which produced a bit more sequence than the other basecallers, and Scrappie events, which produced a bit less. These will explained in the 'Relative read length' section.



### Read identity

<p align="center"><img src="images/read_identity.png" width="90%"></p>

This first analysis tackles the most obvious question: how accurate are the basecalled reads? The plot above shows the read identity distribution, with the median (weighted by read length) marked as a horizontal line. Unaligned reads were given an identity of 0% and fall to the bottom of the distribution. Reads with an actual identity below about 65% usually to fail to align and therefore end up at 0%.

Nanonet performed poorly, with a low median and a significant proportion of unaligned reads. Its curiously high peak of about 99% results from its short output sequences discussed above. While a few Nanonet 'reads' did indeed align to the reference with up to 99% identity, these were really small fragments (hundreds of bp) of larger reads.

Albacore v0.9.1 and Scrappie raw v1.0.0 performed the worst. While their best reads were comparable to other basecallers' best reads, they produced many more reads below 80%. Excluding those versions, Albacore and Scrappie performed well and were comparable to each other. Scrappie raw v1.1.0 rgr_r94 and Albacore v2.0.2 did best and second-best, respectively. Interestingly, Scrappie produced a significant proportion of unalignable reads in each set, whereas Albacore (excluding v0.9.1) did not – I'm not sure why.



### Relative read length

<p align="center"><img src="images/rel_read_length.png" width="90%"></p>

This plot shows the distribution of read length to reference length for each alignment. It shows whether the basecaller is more prone to insertions or deletions. 100% (same length) means that insertions and deletions are equally likely. <100% means that deletions are more common than insertions. >100% means that insertions are more common than deletions.

Curiously, many basecallers produce a distinctly bimodal distribution. For example, Albacore v2.0.2 has a good median value of quite near 100%, but any individual read is more likely to fall either around 99% or 101%. I'm not sure what might be behind this. Albacore v0.9.1 stands out with many overly-long reads, while Scrappie events tends to make short reads. This explains the total yield differences we saw earlier.



### Assembly identity

<p align="center"><img src="images/assembly_identity.png" width="90%"></p>

This analysis is my personal favourite: how accurate are the _consensus_ sequences? I don't particularly care if individual reads have low identity if they can produce an accurate assembly.

Albacore v2.0.2 leads the pack by a decent margin. Most surprisingly, Albacore v0.9.1, which had very poor read identity, does second-best! Scrappie raw v1.0.0 was the worst by far – its assembly was all over the place. I'm not sure what went wrong, but Scrappie raw v1.1.0 seems to have fixed the problem.

<p align="center"><img src="images/rel_assembly_length.png" width="90%"></p>

It's also interesting to look at the assembly relative length, like we did for reads. This shows which basecallers are more prone to consensus sequence insertions (e.g. Albacore v1.0 to v1.2) and which are more prone to deletions (most of the rest).



### Read vs assembly identity

<p align="center"><img src="images/read_assembly_scatter.png" width="90%"></p>

Here I've plotted the median read identity and median assembly identity for all basecallers – zoomed out on the left, zoomed in on the right. The shaded zone is where assembly identity is _worse_ than read identity. That should be impossible (unless you've got a _really_ bad assembler).

This shows how much of each basecaller's error is random vs systematic. If a basecaller had the same read and assembly identities (i.e. on the edge of the shaded zone), that would imply that _all_ of its error is systematic and every read is making the same mistakes. Thankfully, assembly identities are nowhere near that low. Conversely, if a basecaller had an assembly identity of 100%, that would imply a there was very little systematic error so all problems could be ironed out in the consensus.

You might expect that a basecaller's read and assembly identities would be tightly correlated: low-identity reads produce low-identity assemblies and vice versa. That is mostly the case, with Albacore v0.9.1 being the most obvious exception. This suggests that while Albacore v0.9.1 produces reads with a large amount of total error, they have comparatively low systematic error. 



### Post-Nanopolish assembly identity

<p align="center"><img src="images/nanopolish_identity.png" width="90%"></p>

This plot shows the assembly identity distributions after Nanopolish, with pre-Nanopolish assembly identity distributions lightly drawn underneath.

With just two exceptions, all post-Nanopolish assemblies look quite similar. The first exception is Nanonet, where I suspect the truncated reads may have caused problems. The second is Scrappie raw v1.0.0, where the low quality pre-Nanopolish assembly may have been an issue.

The upside seems to be that if you're planning to use Nanopolish, then your basecaller choice may not be very important. Any basecaller, as long as it isn't awful, should be good enough.

The downside is that Nanopolish makes a relatively small improvement to the already good Albacore v2.0.2 assembly. What if a future version of some basecaller can produce a pre-Nanopolish assembly of 99.7% identity or better? I fear that in such a case, Nanopolish may not be able to improve it at all.






## Conclusions

### Recommendations

My current recommendation is simply to use the latest version of Albacore: v2.0.2. It does well on both read and assembly accuracy. Scrappie raw v1.1.0 also did well, especially with the rgr_r94 model. However, since Scrappie is a research product and labelled as a 'technology demonstrator', I think Albacore is a better choice for most users.

The recommendation would have been harder before Albacore v2.0.2 was released. Then, the basecaller with the best assembly accuracy had some of the _worst_ read accuracies: Albacore v0.9.1. Whether or not it would be good choice might depend on your analysis. We've dodged that tough decision for the moment (just use v2.0.2), but may someday be faced with a similar dilemma if a future basecaller excels at consensus accuracy over read accuracy or vice versa.

Finally, Nanonet seems a bit dated and should probably be avoided. However, it does allow for custom training of the neural network and may therefore be of interest to power users interested in experimenting with their own training sets.



### Future work

_My_ future work is easy: trying new versions and new basecallers as they are released and adding them to this analysis. Check back occasionally for new data!

The much harder task lies with the basecaller authors: reducing systematic error. As it currently stands, systematic basecalling errors lead to residual errors in assemblies, even after Nanopolish. This makes it hard to recommend an ONT-only approach for many types of genomics where accuracy matters (read more in [our paper on this topic](http://www.biorxiv.org/content/early/2017/07/07/160614)). If systematic error can be eliminated, then ONT-only assemblies will approach 100% accuracy, and then ONT will be true Illumina alternative.

Did I miss anything important? Can you shed any light on oddities that I couldn't explain? Please let me know through the [issue tracker](https://github.com/rrwick/Basecalling-comparison/issues)!


