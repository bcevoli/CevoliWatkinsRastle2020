# What is semantic diversity and why does it facilitate visual word recognition?
 
This repository contains data, materials and code for computing contextual representations and the semantic diversity metric from Cevoli, Watkins & Rastle (submitted). A pre-print of the manuscript submitted to the journal Behavior Research Methods is available on PsyArXiv.



## Instruction

To replicate the semantic diversity procedure as well as the analysis reported in the pre-print, it is necessary to compute the following steps:
1. Run the *SemanticDiversityProcedure_LSAbasedApproach_Replication.ipynb* notebook found in */procedure* to obtain contextual representations and the semantic diversity metric.
2. Run the *SemanticDiversityProcedure_LSAbasedApproach_Results_Comparison.ipynb* notebook found in */procedure* to compare the results obtained with measures previsouly reported by Hoffman et al. (2013). Please note that in order to produce the figure in this notebook it is necessary to run the code of point 1 multiple times while chnaging the parameters *lemmatization* and *svUnweighting*, which are True and False by deafult respectively. 
3. Run the *A1_ValidationAnalysis.R* script found in */analysis* to replicate the analysis of the newly computed measure of semantic diversity on BLP and ELP.
4. Run the *A1_SimulationAnalysis.R* script found in */analysis* to replicate the analysis of Rodd et al. (2002) and Armstrong & Plaut (2016) on BLP, ELP and semantic diversity. 
5. Run the *SemanticDiversityProcedure_LSAbasedApproach_VisualExploration_AmbiguityCaseStudy.ipynb* notebook found in */procedure* to plot contextual representations of three highly ambiguous words (*calf*, *pupil* and *mole*). 
6. Run the *SemanticDiversityProcedure_LSAbasedApproach_VisualExploration_CorpusMetaData.ipynb* notebook found in */procedure* to plot contextual representations of the entire corpus against metadata.

Eitherwise, it is also possible to run the statistical analyses only as the ouput data of the semantic diveristy procedure are already provided in the */data* folder.




Please note that the first cell each jupyter notebooks contains commented code to help the user install the required packages when necessary (see example below). In order to install a package this cell will need to be uncommented (shortcut to comment/uncomment multiple lines CLRT + / in windows) and the double apostrophes ** after uninstall will need to be replaced with the name of the package. 

```ruby
# Install a pip package in the current Jupyter kernel
import sys
!{sys.executable} -m pip uninstall **
```



## Folder Structure Description 

* /analysis

 This folder contains validation and simulation analysis scripts. 
 The validation analysis script (called *A1_ValidationAnalysis.R*) tests the impact of the newly computed measures of semantic diversity on lexical decision and reading aloud latencies within the English Lexicon Project (ELP; Balota et al., 2007) and the British Lexicon Project (BLP; Keuleers, Lacey, Rastle, & Brysbaert, 2012). The simulation analysis script (called *A2_SimulationAnalysis.R*) replicates Rodd et al. (2002) and Armstrong & Plaut (2016) studies on BLP and ELP lexical decision data, and then simulates the same analysis on the newly computed measure of semantic diversity to test whether this metric can account for lexical ambiguity effetcs in visual word recognition. 

* /data

This folder contains data outputs of the semantic diversity procedure implemented in */procedure* (jupyter notebook named *SemanticDiversityProcedure_LSAbasedApproach_Replication.ipynb*) as well as inputs of the above analysis scripts. 

* /figures

This folder contains figures created in the analysis scripts as well as notebooks comparing semantic diversity results and visually exploring the contextual representations produced by the semantic diversity procedure implemented in */procedure* (jupyter notebook named *SemanticDiversityProcedure_LSAbasedApproach_Replication.ipynb*). These figures can also be found in the manuscript. 

* /procedure
This folder contains jupyter notebooks with the implementation of the Semantic Diveristy Procedure described in Hoffman et al. (2013) and Hsiao & Nation (2018) as well as the semantic diversity results comparison and the visual exploration of the contextual representations. For more information on the semantic diversity implementation see *SemanticDiversityProcedure_LSAbasedApproach_Replication.ipynb*. 

* /stimuli 

This folder contains stimuli from Rodd et al. (2002) and Armstrong & Plaut (2016) with all covariates used in the simulation analysis as well as lexical variables (semantic diversity and frequency) used in the validation analysis. 

* /supplementary

This folder contains summary statistics of the models produced in the analysis scripts. 

