# scRNASeq Pipelines
 contains multiple scripts for data processing and downstream analysis

 1st, download 'scdiff2' analysis software from 
 > 2nd, run 'buildInputExpressionFile5.py' to aggregate 10X genomics output files for pluripotent, cardiac mesoderm, and cardiac progenitor smaple phenotypes and subsequent transformation into '.E' filetypes suitable for further python processing scripts.
 3rd, run prerun pipeline --> 'prerun -i INPUT -o OUTPUT -f FORMAT', where -i represents input scRNA-Seq data that is output from step 2; ...
 4th, run 'anndatapeak.py' to adjust prerun output syntax for final scdiff analysis
 5th, run 'scdiff2' to obtain phenotypic trajectories and evolution analysis.
