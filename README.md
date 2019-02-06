
<p align="center"><img src="images/logo-top.png" alt="logo" width="100%"></p>


# Performance of neural network basecalling tools for Oxford Nanopore sequencing

__Ryan R. Wick<sup>1</sup>, Louise M. Judd<sup>1</sup> and Kathryn E. Holt<sup>1,2</sup>__
<br>
<sub>1. Department of Infectious Diseases, Central Clinical School, Monash University, Melbourne, Victoria 3004, Australia<br>2. London School of Hygiene & Tropical Medicine, London WC1E 7HT, UK</sub>


<p align="center"><img src="images/logo-bottom.png" alt="logo" width="100%"></p>

[![DOI](https://zenodo.org/badge/101518136.svg)](https://zenodo.org/badge/latestdoi/101518136)


This repository contains the scripts used in the preparation of our manuscript on basecalling performance:<br>
[Performance of neural network basecalling tools for Oxford Nanopore sequencing](http://www.biorxiv.org/).

Previous versions of this repository contained the analysis results here in the README, but the current results are now in that manuscript. If you're still interested in the older results, here is a link to the earlier version of this repo: [Comparison of Oxford Nanopore basecalling tools](https://github.com/rrwick/Basecalling-comparison/tree/95bf07476f61cda79e6971f20f48c6ac83e634b3).


## Scripts

[`analysis.sh`](scripts/analysis.sh) is the 'master script' that will run all analyses on a given read set: read-level accuracy, assembly, assembly-level accuracy, nanopolish and nanopolish-level accuracy. It will use the other scripts in its execution. You also might want to edit some of the variables at the start of the script to change things like the output directories and the number of CPU threads. You can also comment out parts of this script if you only want to run some of the analyses.

The `*basecalling*` scripts contain the loops/commands I used to actually run the basecallers. The paths in these scripts are specific to my system – you'll need to edit as appropriate to run them yourself.


## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
