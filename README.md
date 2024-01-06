# scRNASeq Pipelines
 contains multiple scripts for data processing and downstream analysis

 **1st**, download 'scdiff2' analysis software from 'https://github.com/phoenixding/scdiff2/tree/master' <br>
 **2nd**, run 'buildInputExpressionFile5.py' to aggregate 10X genomics output files for pluripotent, cardiac mesoderm, and cardiac progenitor smaple phenotypes and subsequent transformation into '.E' filetypes suitable for further python processing scripts. <br>
 **3rd**, run prerun pipeline --> 'prerun -i INPUT -o OUTPUT -f FORMAT', where '-i' represents input scRNA-Seq data that is output from step 2; '-o' represents the desired output directory to be created that will store subsequently generated data and plots from the run; '-f' represents the data format (raw or normalized). In this case, we will input raw data from the hiPSC-cardiomyocyte differentiation <br> 
 **4th**, run 'anndatapeak.py' to adjust prerun output syntax for final scdiff analysis <br>
 **5th**, run 'scdiff2' pipeline to obtain phenotypic trajectories and evolution analysis --> 'scdiff2.py -i INPUT -o OUTPUT -t TFDNA', where '-i' represents the input scRNA-Seq h5ad data result from the pre-run, '-o' represents the output directory specified in the 4th step, and '-t' represents the transcription factor-DNA interaction data file needed to predict lineage tracing. <br>
