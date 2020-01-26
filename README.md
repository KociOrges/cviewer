#    ```CViewer```
### A Java-based integrated analysis toolkit for whole-genome shotgun metagenomics and other meta’omics technologies

```diff
- Please note that the repository for this project is still under construction!
```

## Introduction to CViewer
The past few years have seen an increased utility of shotgun metagenomics for microbial community surveys over traditional amplicon sequencing. This is made possible by the technological advancement in methods development that enables us now to assemble short sequence reads into longer contiguous regions that can be binned together to identify species they are part of (e.g., through CONCOCT software), and their coding regions can further be annotated against public databases to give an assessment of functional diversity. At the same time, integrated solutions are gaining importance through complementing meta’omics technologies. To consolidate all these realisations on the same sample space, and to fully delineate microbial activity response to environmental factors, it is necessary to include and integrate all levels of gene products, mRNA, protein, metabolites, as well as their interactions in a single platform. Hence, we introduce CViewer, a Java-based software, that integrates with output data from CONCOCT as well as major third party taxonomic and annotation software. The software provides a comprehensive set of multivariate statistical algorithms using the theoretical underpinning of numerical ecology to allow exploratory as well as hypothesis driven analyses. The end result is a highly interactive toolkit with multiple document interface, that makes it easier to unravel useful patterns through Point-and-Click tools whether it is looking at annotated tracks of metagenomic contigs, or exploring enrichments of metabolic pathways and microbial species. 

## Installing
**Step 1: Cloning the repository**

Here you will download a copy of the code and place it somewhere in your directory.
```
$ git clone https://github.com/KociOrges/cviewer.git
```

**Step 2: Unzipping required files**
```
$ cd cviewer
$ unzip Utils.zip
$ unzip Output.zip
```

**Usage**
Once that is done, you can start using CViewer. This can be done by just double-clicking the executable ```CViewer.jar``` file provided in the directory ```cviewer/```.

## Visualisation features

## Clustering, coverage and functional annotation of contigs

CViewer integrates with CONCOCT pipeline (https://github.com/BinPro/CONCOCT) for visualising the contigs. The output from CONCOCT can be loaded into CViewer to inspect the clustering of the contigs across multiple components of the PCA analysis. The software also provides the ability of a focused and interactive exploration of the features of contigs in the context of sample coverages. Additional metadata available on the samples including the treatment groups (say case-control groups) can be uploaded to the software.

The longer sequences of the assembled contigs can further be annotated against public databases to give an assessment of functional diversity. To incorporate this information, CViewer integrates with the output from PROKKA, a tool useful for the genomic annotation of bacterial sequences. The output of the tool is annotation in a GenBank (GBK) file format which serves as an input to the software. One can then use CViewer to explore the annotated genomic features for the contigs of the dataset. For example, one can visualize the CDS regions picked from GBK file, with a means to switch the labels to specific genes or those that were assigned an enzyme identifier. A legend describing summary statistics for the given contig is displayed on the bottom-left part of the panel including length, number of CDS regions, genes and enzyme identifiers.

After the CDS regions are identified using PROKKA, one can search the protein sequences against NCBI’s Conserved Domain Databases (CDD). CDD is a protein annotation resource that consists of a collection of well-annotated multiple sequence alignment models for ancient domains and full-length proteins. Once the CDDs are obtained, we can then use the scripts provided with CViewer to assign them to Clusters of Orthologous Groups (COGs). Each COG consists of a group of proteins found to be orthologous across at least three lineages and likely corresponds to an ancient conserved domain. This annotation process provides an alternative for rapidly describing the functional characteristics of a community of microbes. The COGs identified in the CDS blocks are highlighted by right-clicking and choosing from a pop-up menu. In addition, if there is a need to explore some specific CDS regions further, then a set of them can be selected manually by the user and can be extracted as an output file with their sequence details (say a user wants to analyse how a particular homologous gene differs between different species through third party tools). Finally, integration with statistical utilities for protein sequences allows for the assessment of the structural characteristics of the CDS regions to highlight potential antimicrobial resistant genes.

See animation below:

<img src="https://user-images.githubusercontent.com/30604050/73015452-f5f50e00-3e13-11ea-9110-dd0feb8b86d1.gif" width="900" height="550" />

## KEGG metabolic pathways
PROKKA software also assigns enzyme identifiers to the protein coding regions. The enzyme identifiers are represented by the Enzyme Commission numbers which are assigned by the Nomenclature Committee of the International Union of Biochemistry and Molecular Biology (http://www.chem.qmw.ac.uk/iubmb/enzyme/). An Enzyme Commission number is of the format EC x.y.z.w, where x, y, z, and w are numeric numbers, and describes a numerical classification scheme for the enzymes based on the chemical reactions they catalyse. In CViewer, we have used these numbers to assign the enzymes to molecular level functions by utilizing the Kyoto encyclopaedia of genes and genomes (KEGG) Orthology database. The database consists of a collection of functional orthologs or KEGG Ortholog (KO) groups, that each enzyme is part of, and correspond to the nodes (boxes) of the KEGG pathway maps. The software resolves the KO groups for all the identified enzymes and visualizes KEGG maps with annotation for the mappable enzymes within the pathway. The matched KOs are first converted into presence/absences table for the contigs and combined next with the coverage profiles of that contigs to create microbial pathway abundances across the samples of the metagenomic data.

The complete collection of molecular pathways for metabolism from KEGG database can be visualized and annotated in the software. This includes 11 key metabolic pathways with each of them describing multiple metabolic processes. In each case, the maps are generated on a sample by sample bases i.e., one can visualize a map and inspect the KOs that are upregulated in a sample of choice. Moreover, in addition to providing a visual cue for inspecting the pathways for a given sample, the software allows for comparative analysis of metabolic pathways between samples through downstream statistics, i.e., differential abundance analysis using the Kruskal-Wallis test (described in following section).

<img src="https://user-images.githubusercontent.com/30604050/73015453-f68da480-3e13-11ea-91de-f63ed4552013.gif" width="900" height="550" />

## Phylogenetic diversity
We have considered two  metrics which have been been suggested by Webb (2000) and proved useful for phylogenetic structure analysis, as a measure of phylogenetic clustering or overdispersion, namely the Net Relatedness Index (NRI) and Nearest Taxon Index (NTI). The first focuses on overall clustering, while the latter is more sensitive for capturing the clustering or overdispersion around the tips of the phylogenetic tree. Both metrics become positive as clustering increases and negative with overdispersion. 

To perform phylogenetic analysis with CViewer, one first needs to generate a phylogenetic tree. To do this, the software provides a comprehensive set of essential and highly conserved genes across a wide range of bacteria that can be used for creating a tree, such as those encoding ribosomal proteins, RNA polymerase subunits, and transfer RNA synthetases, with support for manual selection as well. For a more extensive collection of genes useful for phylogenomic analysis the reader is referred to Cambell et al. After the gene has been specified, the gene sequences are extracted for each cluster of the WGS dataset and saved in a .fasta file, along with the average coverages of clusters across the samples of the dataset. The generated sequence .fasta file is next processed with third party tools to generate the phylogenetic tree. Once this is accomplished, the obtained tree and the table describing the abundances of the clusters (generated from the previous step) are uploaded to the software, to calculate the NRI and NTI indices. The user can then choose between two different null models to generate the null communities. These include randomizations within samples while maintaining the species richness or within species while maintaining species frequency. In addition, it is possible to choose from a weighted (where the NRI and NTI metrics for each species are weighted by species abundance) or unweighted approach when performing analysis.


<img src="https://user-images.githubusercontent.com/30604050/73015454-f68da480-3e13-11ea-948d-adb888c553db.gif" width="900" height="550" />


## Alpha diversity
CViewer allows alpha diversity analyses and considers a number of popular indices that can be used to interrogate the input datasets. We have considered the Shannon's (H'), Simpson’s diversity (D1) and its inverse (D2) which account for species richness and abundance. To measure how similar the distributions of species in a community are to each other, the tool provides the Pielou’s evenness. Finally, the relative proportions of the most dominant taxa for a given community dataset can also be explored. See animation below:

![Alpha_diversity 2020-01-12 21_44_37](https://user-images.githubusercontent.com/30604050/72226084-cda02080-3584-11ea-807b-177d882e1d2d.gif)

## Differential abundance
Differential abundance analysis can be valuable when investigating for features (e.g. species/metabolic pathways) that discriminate between multiple treatment groups. CViewer implements the Kruskal-Wallis H statistic for differential analysis and offers a Benjamini-Hochberg correction for multiple comparisons. In addition, pairwise significances between the groups can be explored with the Dunn’s post hoc procedure. The hypotheses to be tested can be uploaded in the software in the form of a .csv file containing hypotheses as headers along with additional metadata on the samples describing the treatment groups. One can then easily navigate through the list of available hypotheses and choose the one to be tested with the Kruskal-Wallis statistic. Finally, significant features are visualized in the form of box plots with appropriate colouring to indicate the treatment groups along with P-values to report significance.

## Correlation
We have considered the Pearson’s product-moment coefficient to measure the degree of the association between two continuous variables, the Kendall’s tau coefficient to compare two quantities measured on at least an ordinal scale and the Spearman’s rank-order correlation as the nonparametric alternative of Pearson’s product-moment coefficient. Prior to analysis, the input tables can be normalized, filtered and/or sorted according to threshold values provided from the user. For a selected correlation coefficient, a full report can be displayed containing the correlation value (R) and the statistical significance (P-value and adjusted P-value) obtained for each pairwise comparison between the elements of the two data matrices.

![Differential_abundance Correlation 2020-01-19 21_09_30](https://user-images.githubusercontent.com/30604050/72688522-0188c200-3b00-11ea-9624-7c862fe4e743.gif)


## Beta diversity
### Principal Component Analysis
Principal component analysis is one of the most popular among the existing dimensionality reduction techniques. The goal here is to find the best summary of the data using a small number of principal components (PCs). We have considered the Nonlinear Iterative Partial Least Squares (NIPALS) for computing PCA, a popular approach in multivariate data analysis that performs well with large datasets.

### Multidimensional Scaling
Similar to PCA, Multidimensional scaling (MDS), also known as Principal Coordinate Analysis (PCoA), is another dimensionality reduction technique that attempts to represent the (dis)similarity between a set of objects in a reduced space, based on their pair-wise dissimilarities given in the form of a distance matrix. The distance matrices can be computed in CViewer choosing from some of the most common distance measures, such as Bray-Curtis, Jaccard and Euclidean distance. For both PCA and MDS, the percentage of variance explained in each dimension is given in the resulting plot.

### Fuzzy Set Ordination
Fuzzy set ordination is used to test effects of pertubation in environmental avariables to community structure. For each of the specified variables, a fuzzy set ordination is calculated and the correlation between the original variable and the fuzzy set is reported. The significance of a particular variable is assessed by comparing a specified threshold p-value and the probability of obtaining a correlation between the data and fuzzy set.

In fuzzy set ordination samples are assigned gradual membership values (fuzzy) ranging from 0 to 1. To achieve this, FSO does not use the raw community data, but rather a similarity (or distance) matrix, which is calculated prior to the statistic. The software supports three different similarity indices to calculate the similarity matrix, namely, the Baroni-Urbani & Buser, Horn and Yule indices. The results are visualised by producing a plot of fuzzy set against original values which is annotated with a correlation between them and a significance label.

### Permutational Multivariate Analysis of Variance
Given a community dataset and a set of predictor physico-chemical variables, PERMANOVA can be used to provide information about the percentage of variation (R2) explained by the given predictors and the significance associated with it (P-value). The results of the PERMANOVA analysis are exported into a summary file containing sources of variation, degrees of freedom, sequential sums of squares, mean squares, K statistics, partial R-squared and P-values, based on N permutations. In addition, the percentage (%) of variation explained (R2) by different predictor variables with annotation for significance can be visualized in the software in the form of a Bar Plot or a Pie chart.

![Beta_diversity 2020-01-12 21_13_52](https://user-images.githubusercontent.com/30604050/72225708-7b5d0080-3580-11ea-9dd0-b2c5809ed7d2.gif)


## Integrated ‘omics analysis
Integrative analysis methods can provide better understanding and interpretation of the composite biological datasets and help elucidate the dynamics of the complex biological systems at a greater scale. CViewer, implements three different techniques that have been proposed for a simultaneous analysis of multiple omics datasets, namely Simultaneous Component Analysis with rotation to COmmon and DIstinctive components (DISCO- SCA), Joint and Individual Variation Explained (JIVE) and Two-way Orthogonal Partial Least Squares (O2PLS). All of them are helpful for providing a “global” view on the biological system of interest by decomposing the variability of the composite omics datasets into a joint variability or common structure, that represents the mechanisms underlying all of the omics datasets under study, and individual variability or distinctive structure, that represents the biological mechanisms underlying a single omics dataset.

### Data pre-processing and model selection
In an integrative or simultaneous analysis of multiple datasets, it is expected that the different blocks comprising the data are linked with a common set of entities. CViewer supports an integrative analysis for two omics datasets that are linked by a common sample space, i.e. the measurements from the different omics data are obtained for the same set of samples (e.g. Crohn’s disease or Healthy individuals).

Prior to analysis, it is useful that the data are pre-processed. As the data are generated from different omics technologies, they may be considerably different in size and describe features that are expressed in scales that are hard to compare, affecting this way the analysis results. For this reason, CViewer provides a number of pre-processing steps that can be applied to correct these differences. When the variables differ largely in scales (or abundance), one may consider centring the variables and/or scale them within each block to a sum of squares of one. In addition, weighting blocks together can be useful to avoid the effects of blocks having considerably different sizes.

Before we choose any of the given methods to perform component analysis, one must first provide the number of components that the composite dataset of omics data is expected to contain, along with their characterization as either common or distinctive. This information is required from each integrative approach (i.e., DISCO-SCA, JIVE and O2PLS) and it is necessary so that the common and individual structures can be successfully identified and separated for each omics dataset. To estimate the number of these common and distinctive components, CViewer provides the Model Selection function. This step is optional, and if this information is known in advance then it can be directly used as an input to the next step of the analysis.

![Multiomics 2020-01-12 23_14_50](https://user-images.githubusercontent.com/30604050/72227185-65a40700-3591-11ea-812e-83802d663581.gif)

## Charts
CViewer implements also a separate tab for the visualisation of charts, with a Line Chart and a Box Plot Chart being supported by the software. This feature is primarily provided as an extra visualisation tool where it is possible to display information as a series of data points, particularly in the case of the Line Chart, and explore trends and patterns between pairs of variables (described from the X and Y axes). For each graph, it is possible to colour the data according to a specific grouping variable of interest. In addition, the Box Plot chart can represent every kind of interval observation using five different measures, i.e. the min, the max, the median and the lower and upper quartiles. This can be a very useful tool when ones need to display the distribution of data, identify outliers, and inspect if and how symmetrical or skewed the data are.

![Charts 2020-01-19 21_33_44](https://user-images.githubusercontent.com/30604050/72688988-87a70780-3b04-11ea-91b7-71c52b2636c4.gif)
