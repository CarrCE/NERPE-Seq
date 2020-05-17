# NERPE-Seq
Non-Enzymatic RNA Primer Extension (NERPE) Deep Sequencing

MATLAB code for planned submission of the manuscript:

Daniel Duzdevich*, Christopher E. Carr* and Jack W. Szostak. Deep sequencing of nonenzymatic RNA primer extension. In preparation. *Joint authors.

## Release notes
Version 1.0 Initial Release

Version 1.1 Updated filter_nerpe_fastq to filter in RNA space

Version 1.2 Updated characterize.m to handle case of no mismatches present

Version 1.3 Updated filter_nerpe_fastq and related code to correct error introduced in version 1.1 where part of the code assumed DNA space. Now all filtering is done in RNA space.

Version 1.4 Added mismatch_context.m and updated characterize.m to handle generating file of context for all mismatches.
