# IPA_pathwayparser
A set of scripts useful for parsing and graphing the results of IPA (Ingenuity Pathway Analysis) exported text files.

This set of scripts provides an easily implemented method for parsing the complete IPA results file, which resembles the following:

```
© 2000-2023 QIAGEN. All rights reserved.

Analysis Details for My Projects->Favorite_project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Workspace	My Projects
Project	Favorite_project
Analysis	deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Date	Mon Jun 19 10:12:58 PDT 2023
Analysis ID	10112233

Analysis Settings
Cutoff before duplicate resolution	true
...
(mol. types = biologic drug OR canonical pathway OR chemical ... OR transmembrane receptor OR transporter)

Canonical Pathways for My Projects->Favorite_project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Ingenuity Canonical Pathways	 -log(p-value)	zScore	Ratio	Molecules
14-3-3-mediated Signaling	-0E00		2.36E-02	PLCG1,PRKCE,TUBB	
Kinetochore Metaphase Signaling Pathway	2.93E01	4.115966043420212	3.33E-01	ANAPC5,AURKB,...,ZWILCH	

Upstream Regulators for My Projects->Favorite_project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Analysis	ID	Upstream Regulator	Expr Log Ratio	Molecule Type	Predicted Activation State	Activation z-score	Notes	Bias Term	Bias-corrected z-score	p-value of overlap	Target molecules in dataset	Mechanistic Network	
deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM	 	CEBPB	0.173	transcription regulator	Activated	9.090	bias	0.441	4.377	1.27E-69	ACTA2,ASF1B,...,XPO1	

Causal Networks for My Projects->Favorite_project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Master Regulator	Molecule Type	Expr Log Ratio	Depth	Predicted Activation State	Notes	Activation z-score	p-value of overlap	Network bias-corrected p-value	...	Relationships Between Master Regulators(Affects/Both)	
ETS1	transcription regulator	-0.224	2	Inhibited	 	-8.560	7.41E-75	1.00E-04	CDKN2A,CSF2,ETS1,FN1	ABCG1,ADAM8,...,ZNF655	142 (4)	4	calcitriol,medroxyprogesterone acetate,NUPR1,...,Vegf	CSF3,FL3,mir-9	 	CASP1,...,ZEB2	CSF1R,IL10,TP53	 	 	 	 	
thimerosal	chemical drug	 	1	 	biased	-1.000	2.70E-02	4.88E-02	thimerosal	BIRC5	1 (1)	1	 	 	 	 	 	 	 	 	 	

Diseases and Bio Functions for My Projects->Favorite_project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Categories	Functions	Diseases or Functions Annotation	p-Value	Predicted Activation State	Activation z-score	Notes	Bias-corrected z-score	Molecules	# Molecules	
Cancer, Hematological Disease, Immunological Disease, Organismal Injury and Abnormalities	MYC mutation positive B-cell lymphoma	MYC mutation positive B-cell lymphoma	9.46E-104	 	 	 	 	TK1,...,CDK1,PLK1	75	

Cancer, Hematological Disease, Immunological Disease, Organismal Injury and Abnormalities	M2 childhood acute myeloid leukemia	M2 childhood acute myeloid leukemia	2.00E-10	 	 	 	 	TOP2A,RRM2,RRM1,IMPDH2,DHFR,POLA1,DCK,POLD1,PRIM1,ANXA1	10	

Tox Functions for My Projects->Favorite_Project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Categories	Functions	Diseases or Functions Annotation	p-Value	Predicted Activation State	Activation z-score	Notes	Bias-corrected z-score	Molecules	# Molecules	
Hepatocellular carcinoma, Liver Hyperplasia/Hyperproliferation	hepatocellular carcinoma	Hepatocellular carcinoma	1.58E-17	Increased	3.228	 	3.491	ABCA13,NEK2,...,CDC25B	142	

Regulator Effects for My Projects->Favorite_Project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Consistency Score	Node Total	Regulator Total	Regulators	Target Total	Target Molecules in Dataset	Disease & Function Total	Diseases & Functions	Known Regulator-Disease/Function Relationships	
4.007	20	1	AREG	18	AURKB,BIRC5,...,SAE1	1	Cell proliferation of tumor cell lines	100% (1/1)	

My Lists for My Projects->Favorite_Project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
My Lists	 -log(p-value)	Molecules 

Tox Lists for My Projects->Favorite_Project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
Ingenuity Toxicity Lists	 -log(p-value)	Molecules 
Cell Cycle: G2/M DNA Damage Checkpoint Regulation	1.24E01	3.08E-01	BRCA1,TOP2A,...,CHEK1

My Pathways for My Projects->Favorite_Project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
My Pathways	 -log(p-value)	Molecules 

Analysis Ready Molecules for My Projects->Favorite_Project->deseq2_TreatmentA_vs_Control_in_Female - 2023-06-19 01:12 PM
ID	Symbol	Expr Intensity/RPKM/FPKM/Counts	Expr Log Ratio	Expr p-value	Expr False Discovery Rate (q-value)
ENSMUSG00000094732.2	1500015L24Rik	6.30516990889	1.39664325030849	2.2629898542166E-8	2.46256155513994E-5

```

### Parsing your IPA results
As a first step, run the IPA parser on your large output file that was exported from IPA.
```
python Parse_IPA.py --help
## Dos filepath
python Parse_IPA.py C:\Users\Molly\ipa\ipa_GroupA_vs_Control_male.txt
## Linux filepath
python Parse_IPA.py /Users/Molly/ipa/ipa_GroupA_vs_Control_male.txt
## Filepaths do not have to be complete, in which case the current directory is used as your output directory.
python Parse_IPA.py ipa_GroupA_vs_Control_male.txt
```
You should then have several output files, as shown below. Note that the default flag "--basic True" can be changed to "--basic False", which would result in additional files being created from the parsed input file.
```
ipa_GroupA_vs_Control_male.txt_canonicalpathways.txt
ipa_GroupA_vs_Control_male.txt_causalnetworks.txt
ipa_GroupA_vs_Control_male.txt_upstreamregulators.txt
```

### Graphing your IPA results
Once the IPA files have been parsed, your can use the content of the provided R script "Automatic_IPA_graphing.R" to plot comparisons between multiple IPA runs. Open your R interface and run the content of the R script interactively. Once all of your input files have been read, you can use the following framework to compare your results:
```
##Designate output directory
out_dir="C:\\Users\\Molly\\Documents\\ipa"
setwd(out_dir)

## Specify your input files
a=list("a.vs.ctl_ko"="C:\\Users\\Molly\\Documents\\ipa\\ipa_treatA.vs.control_knockout.txt",
       "a.vs.ctl_wt"="C:\\Users\\Molly\\Documents\\ipa\\ipa_treatA.vs.control_wildtype.txt",
       "b.vs.ctl_ko"="C:\\Users\\Molly\\Documents\\ipa\\ipa_treatB.vs.control_knockout.txt",
       "b.vs.ctl_wt"="C:\\Users\\Molly\\Documents\\ipa\\ipa_treatB.vs.control_wildtype.txt")
b=lapply(a,function(x) ReadFileForIPAplot(x[1]))

## Specify which IPA results you want to compare
# Example with two datasets
ipa_comparisons_to_make <- c("a.vs.ctl_ko","b.vs.ctl_ko")
gen_plot(ipa_comparisons_to_make,pdf_prefix="Plot_1")
# Example for a four datasets
ipa_comparisons_to_make <- c("a.vs.ctl_ko","a.vs.ctl_wt",
                             "b.vs.ctl_ko","b.vs.ctl_wt")
gen_plot(ipa_comparisons_to_make,k=8,pdf_prefix="Plot_2")
# Example for a single dataset
ipa_comparisons_to_make <- c("a.vs.ctl_ko")
gen_plot(ipa_comparisons_to_make,pdf_prefix="Plot_3")
```

You will now have made some beautiful, if somewhat basic, plots for the following three IPA analysis components: Canonical Pathways, Causal Networks, Upstream Regulators. The current R script is only set up to graph the previously mentioned three analysis components. However, you are welcome to augment the code as needed to accomodate the remaining IPA analysis components (Diseases and Bio Functions, Tox Functions, etc.). Feel free to submit any suggestions/issues for this Github page.



