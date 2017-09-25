# Comparison of Oxford Nanopore basecalling tools

__Ryan R. Wick, Louise M. Judd and Kathryn E. Holt__
<br>
<sub>Department of Biochemistry and Molecular Biology, Bio21 Molecular Science and Biotechnology Institute, University of Melbourne, Australia</sub>


<p align="center"><img src="images/logo.png" alt="logo" width="100%"></p>


## Table of contents

* [Intro](#intro)
* [Data availability](#data-availability)
* [Basecallers tested](#basecallers-tested)
* [Method](#method)
* [Results](#results)
  * [Total yield](#total-yield)
  * [Read identity](#read-identity-1)
  * [Relative read length](#relative-read-length)
  * [Assembly identity](#assembly-identity-1)
  * [Read vs assembly identity](#read-vs-assembly-identity)
  * [Nanopolish assembly identity](#nanopolish-assembly-identity)
* [Conclusions](#conclusions)



## Intro

This repo contains a comparison of available basecallers for Oxford Nanopore Technologies (ONT) sequencing reads. It's a bit like a mini-paper, but I decided to put it here on GitHub (instead of somewhere more papery like [bioRxiv](http://www.biorxiv.org/)) so I can come back and update it as new versions of basecallers are released.

Basecallers, for those not familiar, are the programs which translate the raw electrical signal from an ONT sequencer to a DNA sequence. Basecalling is interesting because it's a hard machine learning problem (modern basecallers all seem to tackle it with neural networks) and because it's a huge part of what makes ONT sequencing good or bad. Getting a piece of DNA through a pore and measuring the current is only half the battle; the other half is in the computer. It's a very active field, with both ONT themselves and independent researchers developing methods.

For each basecaller, I assess the accuracy of the reads and of the resulting assembly. Read accuracy is interesting for obvious reasons – more accurate reads are nice! Assembly accuracy is interesting because shows whether the read errors can 'average out' with depth. In doing so it provides a window into the nature of the basecalling errors. For example, consider a hypothetical set of reads with a mediocre accuracy of 85% but a truly random error profile (i.e. no systematic error). Despite their error rate, these reads could result in a perfect assembly. Now consider a set of reads with an excellent 98% accuracy but they all make the _same mistakes_ (i.e. error is all systematic, not random) – their assembly will also have a 98% error rate. Which read set is better? That probably depends on how you're using them, but in my line of work, I'd prefer the first.

In particular, I hope these results help to answer the question: _Should I go back to old reads and re-basecall them with a newer basecaller?_ Doing so could take a lot of CPU time, so you probably don't want to do it unless it would bring a significant improvement.

As a final note, I used an R9.4 1D dataset of _Klebsiella pneumoniae_ reads for this analysis, so my results may well be biased toward that kind of data. I'm not sure how consistent these results are with other data types, e.g. eukaryote genomes, R9.5 flow cells, 1D<sup>2</sup> kits, etc.





## Data availability

If you'd like to try this analysis using the same data, here are the relevant links:
* [Reference hybrid assembly](https://figshare.com/articles/Unicycler_v0_4_0_assemblies_hybrid_Illumina_and_ONT_/5170750) (barcode01.fasta.gz)
* [Raw fast5 files](https://figshare.com/articles/Raw_ONT_reads_-_barcode_1/5353210)
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



### Chiron

[Chiron](https://github.com/haotianteng/chiron) is a third-party basecaller developed by [Haotian Teng](https://github.com/haotianteng) and others in [Lachlan Coin's group](https://imb.uq.edu.au/profile/647/lachlan-coin) at the University of Queensland.

```
chiron call -i raw_fast5 -o chiron_v0.2 --batch_size 1000
```

While testing Chiron, I noticed a curious effect. The `--batch_size` parameter was described as controlling performance: a larger value improves performance but increases RAM requirements. While I found this to be true, I discovered a secondary effect: larger `--batch_size` values also improved the read accuracy. The default value is `100` and I saw modest read accuracy improvements up to about `500`, after which accuracy plateaued. In my tests I used `--batch_size 1000`for both performance and accuracy.

NOTE: I'm still running Chiron tests, so they aren't included (yet) in the results below. Check back soon!



### basecRAWller

[basecRAWller](https://basecrawller.lbl.gov/) is a third-party basecaller developed by Marcus Stoiber and James Brown at the Lawrence Berkeley National Laboratory. Unlike other basecallers, it focuses on _streaming_ basecalling. I.e. it can basecall data using any part of a read and does not require the entire sequence to be present. This could potentially be used to basecall reads even before they are finished sequencing. As discussed in the [basecRAWller paper](https://www.biorxiv.org/content/early/2017/05/01/133058), the ability to perform streaming basecalling does negatively impact accuracy.

```
basecRAWller call --fast5-basedirs raw_fast5 --out-filename basecrawller_v0.1.fasta
```

NOTE: I'm still running basecRAWller tests, so they aren't included (yet) in the results below. Check back soon!



### Excluded basecallers

I am also interested in [Guppy](https://github.com/nanoporetech/guppy), an ONT basecaller designed for fast GPU-accelerated performance. However, it is only available to users who have signed the ONT developer license agreement (which is why you might have gotten a 404 if you just tried that link). If and when Guppy becomes publicly available, I'll add it to this comparison too.

Unfortunately, I cannot compare with the old cloud-based Metrichor basecalling, as it's no longer available. I also cannot test the integrated basecalling in MinKNOW (ONT's sequencing software). I believe MinKNOW's integrated basecalling shares much in common with Albacore, but I don't know which Albacore versions correspond to which MinKNOW versions.





## Method

If you'd like to try this analysis for yourself, here's what you need to do. I tried to make this process somewhat flexible, but some aspects may be particular to my setup, so you'll probably need to modify some parts to make it work for you. In particular, the `nanopolish_slurm_wrapper.py` script assumes you're using a SLURM-managed cluster, so other users will probably need to change that one.



### Required files

You'll obviously need a set of ONT reads. Put them in a directory named `01_raw_fast5`. I used the same ones from our recent paper: [Completing bacterial genome assemblies with multiplex MinION
sequencing](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000132). Check out that paper if you're in the wet lab side of things.

You'll also need Illumina reads for the sample (`illumina_1.fastq.gz` and `illumina_2.fastq.gz`) and a good reference sequence (`reference.fasta`), e.g. a completed hybrid assembly.

My reads came from a barcoded run, so I first had to collect only the fast5 files for my sample. I did this by analysing the fastq file of our confidently-binned reads (again, see [the paper](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000132) for more info). This process should have excluded most of the very low quality reads, because such reads would not have been confidently binned. I also tossed out any fast5 files less than 100 kilobytes in size to remove shorter reads, though this step may not be necessary.



### Required tools

The following tools must be installed and available in your `PATH`: [minimap2](https://github.com/lh3/minimap2), [Filtlong](https://github.com/rrwick/Filtlong), [Porechop](https://github.com/rrwick/Porechop), [Racon](https://github.com/isovic/racon), [Rebaler](https://github.com/rrwick/Rebaler), Nanopolish and SAMtools.

I used [Nanopolish](https://github.com/jts/nanopolish) [v0.8.1](https://github.com/jts/nanopolish/releases/tag/v0.8.1) to improve each assembly and then assessed identity again. Since v0.8.0, Nanopolish can be run without event-data-containing fast5 files, which lets me run it with any basecaller! However, for non-Albacore basecallers I did have to alter read names so they were Albacore-like and compatible with the `nanopolish index` command.



### Basecalling

The specifics here depends on the basecaller – the commands I used are described above. When basecalled reads are ready, put them in a `02_basecalled_reads` directory as either `*.fastq.gz` or `*.fasta.gz` files.



### Read ID to fast5 file

It is also necessary to make a `read_id_to_fast5` file which contains two columns: read_ID in the first and fast5 filename in the second. This is to ensure that all basecalled read files are named consistently and in a format compatible with Nanopolish.

For example:
```
0000974e-e5b3-4fc2-8fa5-af721637e66c_Basecall_1D_template	5210_N125509_20170425_FN2002039725_MN19691_sequencing_run_klebs_033_restart_87298_ch173_read25236_strand.fast5
00019174-2937-4e85-b104-0e524d8a7ba7_Basecall_1D_template	5210_N125509_20170424_FN2002039725_MN19691_sequencing_run_klebs_033_75349_ch85_read2360_strand.fast5
000196f6-6041-49a5-9724-77e9d117edbe_Basecall_1D_template	5210_N125509_20170425_FN2002039725_MN19691_sequencing_run_klebs_033_restart_87298_ch200_read1975_strand.fast5
```



### Run analysis

The `analysis.sh` script automates most of the remaining steps. It will:
1) Change the read names to a consistent, Nanopolish-friendly format (`fix_read_names.py`).
2) Align the reads to a reference and make a tsv file of read accuracies (`read_length_identity.py`). This only uses the aligned parts of the read to calculate the read's identity. The definition used for 'identity' is the same as how BLAST defines it: the number of matching bases in the alignment divided by the total bases in the alignment (including gaps). If less than 50% of the read aligned, it is deemed unaligned and given an identity of 0%. This script also determines the read length to reference length ratio for each read, to see if insertions or deletions are more likely.
3) Prepare reads for assembly (`porechop` and `filtlong`).
4) Do a reference-based assembly (`rebaler`).  [Rebaler](https://github.com/rrwick/Rebaler), conducts multiple rounds of [Racon](https://github.com/isovic/racon), so the final assembly accuracy is defined by the Racon consensus (which in my experience is a bit higher accuracy than a [Canu](https://github.com/marbl/canu) assembly).
5) Assess the assembly in the same manner as it did for reads (`chop_up_assembly.py` and `read_length_identity.py`). By chopping the assembly into 10 kbp pieces and assessing them like reads, we can get a _distribution_ of assembly identity instead of just a single value.
6) Run Nanopolish (`nanopolish_slurm_wrapper.py`).
7) Assess the Nanopolished assembly (`chop_up_assembly.py` and `read_length_identity.py`).



### Generate figures

Put all of your resulting tsv files in a `results` directory and run `plot_results.R` to generate figures.





## Results

### Total yield

<p align="center"><img src="images/total_yield.png" width="90%"></p>

You might expect that each basecaller would produce approximately the same total yield. E.g. a read that makes a 10 kbp sequence in one basecaller would be about 10 kbp in each basecaller. That's true for the most part, but Nanonet is a notable exception. For most reads, it produced a much shorter sequence than other basecallers, sometimes drastically so. For example, all versions of Albacore basecalled one read (d2e65643-98b4-4a61-ad36-119e55250b28) to a 34+ kbp sequence. Nanonet produced 518 bp for the same read. I don't have an explanation for this odd behaviour.

Other oddities you might notice are Albacore v0.9.1, which produced a bit more sequence than the other basecallers, and Scrappie events, which produced a bit less. These will explained in the 'Relative read length' section.



### Read identity

<p align="center"><img src="images/read_identity.png" width="90%"></p>

This first analysis tackles the most obvious question: how accurate are the basecalled reads? The plot above shows the read identity distribution, with the median (weighted by read length) marked as a horizontal line. Unaligned reads were given an identity of 0% and fall to the bottom of the distribution. Reads with an actual identity below about 65% usually to fail to align and therefore end up at 0%.

Nanonet performed poorly, with a low median and a significant proportion of unaligned reads. Its curiously high peak of about 99% results from its short output sequences discussed above. While a few Nanonet 'reads' did indeed align to the reference with up to 99% identity, these were actually just small fragments (hundreds of bp) of larger reads.

Albacore v0.9.1 and Scrappie raw v1.0.0 had the lowest median identities. While their best reads were comparable to other basecallers' best reads, they produced many more reads below 80%. Excluding those versions, Albacore and Scrappie performed well and were comparable to each other. Scrappie raw v1.1.0 rgr_r94 and Albacore v2.0.2 did best and second-best, respectively. Interestingly, Scrappie produced a significant proportion of unalignable reads in each set, whereas Albacore (excluding v0.9.1) did not – I'm not sure why.



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



### Nanopolish assembly identity

<p align="center"><img src="images/nanopolish_identity.png" width="90%"></p>

This plot shows the assembly identity distributions after Nanopolish, with pre-Nanopolish assembly identity distributions lightly drawn underneath.

In every case, Nanopolish improved the assembly accuracy, and with just two exceptions, all post-Nanopolish assemblies look quite similar. The first exception is Nanonet, where I suspect the truncated reads may have caused problems. The second is Scrappie raw v1.0.0, where the low quality pre-Nanopolish assembly may have been an issue.

The upside seems to be that if you're planning to use Nanopolish, then your basecaller choice may not be very important. Any basecaller, as long as it isn't awful, should be fine. The downside is that Nanopolish makes a relatively small improvement to the already good Albacore v2.0.2 assembly. What if a future basecaller can produce a pre-Nanopolish assembly of 99.7% identity or better? I fear that in such a case, Nanopolish may not be able to improve it at all.





## Conclusions

### Recommendations

My current recommendation is simply to use the latest version of Albacore: v2.0.2. It does well on both read and assembly accuracy. Scrappie raw v1.1.0 also did well, especially with the rgr_r94 model. However, since Scrappie is a research product and labelled as a 'technology demonstrator', I think Albacore is a better choice for most users.

The recommendation would have been harder before Albacore v2.0.2 was released. Then, the basecaller with the best assembly accuracy had some of the _worst_ read accuracies: Albacore v0.9.1. Whether or not it would be good choice might depend on your analysis. We've dodged that tough decision for the moment (just use v2.0.2), but may someday be faced with a similar dilemma if a future basecaller excels at consensus accuracy over read accuracy or vice versa.

Finally, Nanonet seems a bit dated and should probably be avoided. However, it does allow for custom training of the neural network and may therefore be of interest to power users interested in experimenting with their own training sets. That being said, ONT has open-sourced their [Sloika project](https://github.com/nanoporetech/sloika) for training Albacore/Scrappie neural nets, so users doing custom training may want to use that over Nanonet.



### Training and methylation

The training set used for any given basecaller may have a huge impact on the quality of the results. [Tim Massingham](https://github.com/tmassingham-ont) (the author of Scrappie) mentioned [here](https://github.com/rrwick/Basecalling-comparison/issues/1) that Albacore v2.0.2 and Scrappie raw v1.1.0 rgrgr_r94 are very similar. The key difference is that Scrappie was trained with a human-only dataset, which possibly explains why Albacore v2.0.2 did much better for my _Klebsiella_ sample.

Methylation may also be an important factor related to training, as a methylated base would be expected to produce a different signal in the pore than its unmethylated counterpart. If a basecaller was trained with mostly or entirely unmethylated DNA, then we might expect it to give low accuracy on heavily methylated samples. This may explain much of the residual error I saw in my assemblies. If this is the case, two solutions jump to mind: 1) PCR your DNA before ONT sequencing, or 2) train the basecaller's neural network using DNA (_with_ methylation) which is more like your sample. I don't like the first solution (more wet lab work), but the second seems promising. If it really is as simple as using a better training set, then significant basecaller improvements are potentially just around the corner.



### Future work

_My_ future work is easy: trying new versions and new basecallers as they are released and adding them to this analysis. Check back occasionally for new data!

The much harder task lies with the basecaller authors: reducing systematic error. As it currently stands, systematic basecalling errors lead to residual errors in assemblies, even after Nanopolish. This makes it hard to recommend an ONT-only approach for many types of genomics where accuracy matters (read more in [our paper on this topic](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000132)). If systematic error can be eliminated, ONT-only assemblies will approach 100% accuracy, and then ONT will be a true Illumina alternative.

Did I miss anything important? Can you shed any light on oddities that I couldn't explain? Please let me know through the [issue tracker](https://github.com/rrwick/Basecalling-comparison/issues)!
