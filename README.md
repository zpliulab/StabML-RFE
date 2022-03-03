# [StabML-RFE (stable machine learning-recursive feature elimination method)](https://github.com/zpliulab/StabML-RFE)

![Screenshot](Data/HGSOC.jpg)

In this work, we provide **a stable machine learning-recursive feature elimination method** named **StabML-RFE ** for identifying **robust diagnostic biomarkers** of high-grade serous ovarian cancer (HGSOC) based on gene expression data. The successful identification of HGSOC biomarkers will be beneficial to reduce the risk of ovarian cancer among women for early disease detection. **Obviously, the proposed method of discovering biomarkers for HGSOC can be easily extended for other complex diseases**.


## StabML-RFE
<!--START_SECTION:news-->
* **StabML-RFE **: A method of stable machine learning-recursive feature elimination for identifying robust biomarker from gene expression data. 
* If you have any questions about **StabML-RFE **, please directly contact the corresponding author [Prof. Zhi-Ping Liu](https://scholar.google.com/citations?user=zkBXb_kAAAAJ&hl=zh-CN&oi=ao) with the E-mail: zpliu@sdu.edu.cn
<!--END_SECTION:news-->


## Citation
Li, Lingyu, and Zhi-Ping Liu. "**Identifying robust diagnostic biomarkers of high-grade serous ovarian cancer by stable machine learning-recursive feature elimination method**." Submit to Journal. 


## Data
<!--START_SECTION:news-->
* In the **Datanew**, **DataPythonTop** and **DataPythonTop20** files, we only give some necessary files by each R or Python program. 
* Some of these input files only give the first few lines, but this does not affect the results of the work (**StabML-RFE**).
* File **Datanew**, it contains most files that the codes need to input or the codes output.
* Files **DataPythonTop** and **DataPythonTop20**, they respectively contains the eight ranking lists obtained from eight ML-RFE methods and the optimal feature sets induced by the top-ranked 20 features of each method.
<!--END_SECTION:news-->


## R code for StabFS 
The **serial number (1) (2) ... (10)** represents the order in which the program runs in our work.
<!--START_SECTION:news-->
* (1) GSE69284_expr.R ---- Processing original data of GSE69284 and get the (scale) data with corresponding labels.
* (2) DE_gene.R ---- Identify the differentially expressed genes (DEGs) on a dataset, find candidates with adjusted adj.P.Val < 0.05 and |logFC|>1.5, and get the DEGs’s expression data file “matrix_DE” used to input Python codes to conduct ML-RFE progression.
* (3) rfeExample.py ---- Perform the seven ML-RFE methods with the random seed p = 2022, which is sued to keep the same result in each run. Thus the ranking files of all DE features are obtain (ranklist_SVMrfe.txt, ranklist_ABrfe, ranklist_RFrfe.txt, ranklist_DTrfe.txt, ranklist_GBDTrfe.txt, ranklist_XGBrfe.txt, ranklist_NBrfe.txt). 
* (4) rfe_nnet.R ---- Conduct NNET-RFE method and obtain “ranklist_NNETrfe.txt”.
* (5) RFE_ROC_newTest.R ---- Seleted the top-ranked 20 features of each list obtained from eight ML-RFE methods, take them as the optimal feature subsets, and test the classification results of each optimal feature subsets.
* (6) stability_HGSOCpythonTop.R ---- Calculate the stability index of all combinations for quantifying the similarity of two or more feature sets, and identify the potential biomarkers. It proves the stability/robustness of the identified biomarkers selected by our subset identification strategy of different machine learning and feature selection methods.
* (7) feature_study.R -- Extract the expression value of the biomarkers from the discovery dataset. 
* (8) cornew.R  --  It mainly make internal validation including the differential analysis (Plot the differential barplot and label the significance of the identified biomarkers on the discovery dataset), heatmap plot, the correlation of the biomarkers and the up/down regulation, and the network of 18 biomarkers with correlation larger than 0.8.
* (9) cluster_HGSOC.R  --  Make GO and pathway functional enrichment analysis of the identified biomarkers.
* (10) feature_selectnew.R ---- Extract data from independent dataset and origina discovery dataset based on identified biomarker.
* (11) class_featurenew.R " ---- Verify the identified biomarkers on 5 independent dataset and obtain the true label and predict label.
* (12) ROC_plotnew.R" ---- Plot ROC curves and AUC value of 5 independent validation datasets  and plot ROC curves.
* (13) DEplot_add.R ---- Show the classification performances of the external validation in the form of bar plots.
<!--END_SECTION:news-->


## Stability of Feature Selection Techniques for Bioinformatics Data
<!--START_SECTION:news-->
* [Feature selection is one of the most fundamental problems in data analysis, machine learning, and data mining](https://doi.org/10.1007/978-3-030-64583-0_19). Especially in domains where the chosen features are subject to further experimental research, the stability of the feature selection is very important. Stable feature selection means that the set of selected features is robust with respect to different data sets from the same data generating distribution.
* For data sets with similar features, the evaluation of feature selection stability is more difficult. An example of such data sets is gene expression data sets, where genes of the same biological processes are often highly positively correlated.  Here, stability measures that take into account the similarities between feature subsets are defined as **Hamming stability**.
<!--END_SECTION:news-->


## LogReg (2022), Zhi-Ping Liu all rights reserved
This program package is supported by the copyright owners and coders "as is" and without warranty of any kind, express or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. In no event shall the copyright owner or contributor be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, without limitation, procurement of substitute goods or services; loss of use, data, or profits; or business interruption), regardless of the theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) for any use of the software, even if advised of the possibility of such damages.
