# NERPE-Seq
Non-Enzymatic RNA Primer Extension (NERPE) Deep Sequencing

Citing NERPE-Seq: Please use:

Daniel Duzdevich*, Christopher E. Carr* and Jack W. Szostak. Deep sequencing of nonenzymatic RNA primer extension. Nucleic Acids Research, gkaa400, https://doi.org/10.1093/nar/gkaa400, 19 May 2020. Preprint: bioRxiv 10.1101/2020.02.18.955120 *Joint authors. Code: https://github.com/CarrCE/NERPE-Seq 

Additional follow-on work has been published as:

Daniel Duzdevich, Christopher E Carr, Dian Ding, Stephanie J Zhang, Travis S Walton, Jack W Szostak, Competition between bridged dinucleotides and activated mononucleotides determines the error frequency of nonenzymatic RNA primer extension, Nucleic Acids Research, Volume 49, Issue 7, 19 April 2021, Pages 3681–3691, https://doi.org/10.1093/nar/gkab173

## Release notes
Version 1.0 Initial Release

Version 1.1 Updated filter_nerpe_fastq to filter in RNA space

Version 1.2 Updated characterize.m to handle case of no mismatches present

Version 1.3 Updated filter_nerpe_fastq and related code to correct error introduced in version 1.1 where part of the code assumed DNA space. Now all filtering is done in RNA space.

Version 1.4 Added mismatch_context.m and updated characterize.m to handle generating file of context for all mismatches.

Version 1.5 Updated characterize.m to also generate transition map of complementary set in addition to all products/templates; updated transition_map.m to improve transition map figure with nulls removed.

Version 1.6 Updated seqspace_cube.m to use rendering for non-pixelated output.

Version 1.7-1.9 Updated figures to auto-close to facilitate running many samples.

Version 2.0 Illumina sequencing read length was restricting ability to identify products at 13 bases or longer. The result was that many fix2 and fix3 sequences were truncated in the Illumina reads. To adjust for this, this updated ProcessSamples allows additional columns to be added to the input excel file, e.g. columns "fix2" and "fix3" to define shorter fix2 and fix3 sequences. In this way, longer products can be identified up to 18 bases long when using the following:

fix2 = 'AGATCGGAAGAGCACAC'
fix3 = 'GATCGTCGGACTGTAGA'

instead of the defaults:

fix2 = 'AGATCGGAAGAGCACACGTCTGA'
fix3 = 'GATCGTCGGACTGTAGAACTCTG'

