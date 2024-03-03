<img width="280" height="130" alt="Screenshot 2020-08-10 at 11 30 22" src="https://user-images.githubusercontent.com/30604050/90171121-eed62300-dd98-11ea-9ef2-0336eed842b2.png">

### A Java-based statistical framework for integration of shotgun meta-genomics with other omics datasets


#### CViewer is developed by Orges Koci as part of his PhD under the supervision of Dr Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz/) in his Environmental Omics Lab ([Water & Environment Group](https://www.gla.ac.uk/schools/engineering/research/divisions/i_e/researchthemes/w_e/)) and jointly with the [BINGO Group](https://www.gla.ac.uk/schools/medicine/research/humannutrition/thebingogroup/) with the critical feedbacks by Profs Konstantinos Gerasimidis and Richard Russell who served in his thesis committee. 

```diff
- Please note that this is a beta version and the repository for this project is still under construction.
```

## Publication
- A preprint of this work is now available on bioRxiv: [https://www.biorxiv.org/content/10.1101/2023.06.07.544017v1](https://www.biorxiv.org/content/10.1101/2023.06.07.544017v1)

## Table of Contents
- [Introduction to CViewer](#introduction-to-cviewer)
- [Dependencies](#dependencies)
- [Installation & Usage](#installation--usage)
  - [For Intel-based Users (Standard Usage)](#for-intel-based-users-standard-usage)
  - [For ARM-based (Apple Silicon) Macs](#for-arm-based-apple-silicon-macs)
- [Visualizing the Metagenomics Contigs](#visualizing-the-metagenomics-contigs)
- [Functional Annotation of the Metagenomics Contigs](#functional-annotation-of-the-metagenomics-contigs)
- [KEGG Metabolic Pathways](#kegg-metabolic-pathways)
- [Phylogenetic Diversity](#phylogenetic-diversity)
- [Data Analyses](#data-analyses)
  - [Alpha Diversity](#alpha-diversity)
  - [Beta Diversity](#beta-diversity)
    - [Principal Component Analysis](#principal-component-analysis)
    - [Multidimensional Scaling](#multidimensional-scaling)
    - [Fuzzy Set Ordination](#fuzzy-set-ordination)
    - [Permutational Multivariate Analysis of Variance](#permutational-multivariate-analysis-of-variance)
- [Differential Abundance](#differential-abundance)
  - [Kruskal-Wallis H](#kruskal-wallis-h)
  - [Friedman Test](#friedman-test)
- [Correlation](#correlation)
- [Integrated ‘omics Analysis](#integrated-omics-analysis)
- [Charts](#charts)
- [References](#references)


## Introduction to CViewer
The past few years have seen an increased utility of shotgun metagenomics for microbial community surveys over traditional amplicon sequencing. This is made possible by the technological advancement in methods development that enables us now to assemble short sequence reads into longer contiguous regions that can be binned together to identify species they are part of, e.g., through CONCOCT[[1]](#1). The advantage of shotgun metagenomics is that  coding regions of these contigs can further be annotated against public databases to give an assessment of the functional diversity. With integrated solutions gaining importance by complementing metagenomics with other meta’omics technologies (e.g., metabolomics), there is a need to have a single platform to consolidate all these realisations on the same sample space. Thus, we have developed ```CViewer``` which aims to integrate all levels of gene products, mRNA, protein, metabolites for microbial communities and allows exploration of their response to environmental factors through multivariate statistical analysis. 
 
The current version of ```CViewer``` is a highly interactive Java application tailored to the needs of non-expert users. It integrates with the output from CONCOCT as well as major third party taxonomic (KRAKEN[[2]](#2) etc.) and annotation software (PROKKA[[3]](#3)). The theoretical underpinning of numerical ecology allows exploratory as well as hypothesis driven analyses, emphasizing on functional traits of microbial communities and phylogenetic‐based approaches to community assembly, particularly abiotic filtering. 

## Dependencies
- ```Operating system:``` Platform independent
- ```Requirements:``` Java Platform Standard Edition 1.8.0 or later
- ```License:``` MIT

## Installation & usage
**Step 1: Cloning the repository**

Here, you will download a copy of the code and place it somewhere in your directory.
```
$ git clone https://github.com/KociOrges/cviewer.git
```

**Step 2: Unzipping required files**
```
$ cd cviewer
$ unzip Utils.zip
$ unzip Output.zip
```

### For Intel-based Users (Standard Usage)

Users with Intel-based processors can start CViewer by simply double-clicking the ```CViewer.jar``` file located in the ```cviewer/``` directory. Alternatively, for those who prefer using the terminal, execute the following command within the ```cviewer/``` directory:

```
java -jar CViewer.jar
```

### For ARM-based (Apple Silicon) Macs

Apple Silicon Mac users (e.g., M1 chip) need to take a few extra steps due to architectural differences. This involves unzipping the dependencies.zip file and ensuring you have the compatible JavaFX SDK and Java version for your machine.

****1. Unzipping Additional Dependencies****

Navigate to the cviewer/ directory and unzip the **dependencies.zip** file to include all necessary JavaFX modules and libraries:
```
unzip dependencies.zip
```

****2. Running the Software****

To execute CViewer on Apple Silicon Macs, use the following command, which specifies the necessary JavaFX module path and includes all dependencies:

```
java --module-path "./dependencies/javafx-sdk-21.0.1" --add-modules javafx.controls,javafx.fxml,org.controlsfx.controls,javafx.swing -cp "CViewer.jar:./dependencies/*:./dependencies/." MainGUI_FX
```

These dependencies are provided as a workaround for convenience and were tested on an Apple M1 Max with Java version "21.0.1" and JavaFX SDK 21.0.1.

****3. Customizing for Your System****

We recommend downloading the JavaFX SDK and Java version that best match your system's requirements. Ensure the JavaFX version is compatible with your Java version. You can download the JavaFX SDK from [Gluon](https://gluonhq.com/products/javafx/).

Should you need a different version of JavaFX SDK, download it and place it within the dependencies folder. Update the execution command to reflect the path to your downloaded JavaFX SDK version.

****4. Note on Security Permissions****

Depending on your system's security settings, you may need to adjust permissions to allow execution of the CViewer.jar file. If prompted, confirm that you trust the source of this file.


## Visualizing the metagenomics contigs

```CViewer``` integrates with CONCOCT pipeline (https://github.com/BinPro/CONCOCT) for visualising the contigs. The output from CONCOCT can be loaded into the software to inspect the clustering of the contigs across multiple components of the PCA analysis in the context of samples coverages. To explore the contigs in CViewer, first, we will need to obtain the required files from CONCOCT and then import it into the software. Here, as an example, we have put the generated files inside the folder *example_datasets/*. 

**Step 1: Collecting the required data**

For this tutorial, you will need the files listed below:
* ```PCA_transformed_data_gt1000.csv:``` the *N x D* CSV file that contains the *D* PCA components for each contig 
* ```concoct_inputtable.csv:``` the *N x S* TSV file that contains the information for the coverage values for the *S* samples of the dataset
* ```clustering_gt1000.csv:``` the *N x 1* CSV file that gives the labelling of which cluster each contig belongs to

In addition, we will need the file with the taxonomic labels for the clusters of contigs which, in this example, is obtained from Kraken software.

* ```taxonomy.csv```


**Step 2: Importing the data into CViewer**

To start importing the above information into CViewer, you will need to navigate youself to the "Main" tab. Next, look for the "Input" section and click on the arrow on the left-hand side to expand (or collapse) if necessary, the section's menu. CViewer provides a set of self-explanatory buttons for importing the required data. Note that each button is provided for importing a specific file at a time. This is described below:

Click the **Open** button next to
* ```PCA:``` to import the file that contains the PCA components (e.g. PCA_transformed_data_gt1000.csv)
* ```Coverage:``` for the file that contains the information for the coverage values (e.g. concoct_inputtable.csv)
* ```Clustering:``` for the file that gives the labelling of which cluster each contig belongs to (e.g. clustering_gt1000.csv)
* ```Taxonomy:``` for the file that contains the taxonomic information for the dataset (e.g. taxonomy.csv)

This is also illustrated in the following animation: [![animation](https://user-images.githubusercontent.com/30604050/88484527-06bf4180-cf67-11ea-9388-a9890a91dc13.png)](https://www.youtube.com/watch?v=DmD4qEQ71mc&feature=youtu.be)


**Step 3: Exploring the contigs through PCA plots**

Once we have imported the above data into the software we can start exploring the metagenomics contigs. To do so, you will need to click and expand the ```Principal Component Analysis``` section and press the ```Apply``` button under the ```Components``` _slider-bar_. This will initiate the processing of the input files and within a few minutes you should see a PCA plot on your screen. It is possible then to inspect the clustering of the contigs further up to eight dimensions of the PCA analysis while one can also explore the features of contigs in the context of sample coverages. The user can choose to switch the view of the PCA panel to a ‘coverage’ mode and inspect the coverage information for the contigs across the samples of the dataset. The size of the contigs in the PCA space shows the coverages of the contigs and serves as a visual cue for the species that are abundant in a given sample with the _slider-bar_ giving a means to shift between samples. Additional metadata associated with the samples can be imported from a file, in order to explore different treatment groups (e.g., samples for Crohn's Disease vs. healthy individuals). Set of buttons provides additional functionalities by showing/hiding the legends for the dataset clusters (cluster number, taxonomy etc.) or adjusting the font size. More particularly, the software provides the following options:

* ```Show legend:``` whether legend with information on the dataset clusters (cluster number, taxonomy etc.) should be displayed
* ```Show taxonomy:``` whether the taxonomic labels for the dataset clusters should be displayed (if applicable)
* ***```+```***, ***```-```*** buttons: used to increase/ decrease legend's font size
* ```Show coverage:``` whether ‘coverage’ mode should be activated to display contigs coverages
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

* ```Components:``` slider-bar showing the dimensions of the PCA analysis, used for selecting and inspecting the variations of contigs clusters across the available PCA components (enabled when "Show coverage"=false). Window is updated by pressing the ```Apply``` button.
* ```Samples``` button: used for importing additional metadata associated with the samples from an external file and for selecting the labels for the dataset samples
* ```Samples``` slider-bar (under the ```Samples``` button): used for selecting the sample for which we wish to inspect the contigs coverages (enabled when "Show coverage"=true). Window is updated by pressing the ```Apply``` button.

The above functionalities are demonstrated in the video right below:
[![animation](https://user-images.githubusercontent.com/30604050/88486314-c36bcf80-cf74-11ea-8d75-4955cfa32a8d.png)](https://www.youtube.com/watch?v=nNBwd6isETQ&feature=youtu.be)

## Functional annotation of the metagenomics contigs

One can explore the annotated genomic features for the contigs of the dataset by visualising the CDS regions picked from GenBank file (given from PROKKA which, is already run on the contigs in a prior step, https://github.com/tseemann/prokka), with a means to switch the labels to specific genes or those that were assigned an enzyme identifier. One can also search the identified protein sequences against NCBI’s **[Conserved Domain Databases (CDD)](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml)**. CDD is a protein annotation resource that consists of a collection of well-annotated multiple sequence alignment models for ancient domains and full-length proteins. Once the CDDs are obtained, we can then use the scripts provided with ```CViewer``` to assign them to **Clusters of Orthologous Groups (COGs)**. Each COG consists of a group of proteins found to be orthologous across at least three lineages and likely corresponds to an ancient conserved domain. This annotation process provides an alternative for rapidly describing the functional characteristics of a community of microbes. To explore this information in CViewer, we first need to import the required data obtained from PROKKA and NCBI's CDD database located inside the folder *example_datasets/* and use the "Input" section of the "Main" tab again in the similar manner as previously described:

**Step 1: Collecting the data**

For this part of the analysis, you will need the following files:
* ```PROKKA_XXXXXXX.gbk:``` the file that contains all the annotation tracks for the *N* contigs of the dataset, obtained from PROKKA software (e.g. example_datasets/PROKKA.gbk)
* ```PROKKA_XXXXXXX.faa:``` the file describing all of the predicted protein coding sequences for the *N* contigs of the dataset (e.g. example_datasets/PROKKA.faa)

To obtain the COGs, we will first use the [RPS-BLAST](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rpsblast/) to find COG annotation for the protein sequences using the script RPSBLAST.sh that is provided with the CONCOCT software. A complete example for doing this is given [here](https://concoct.readthedocs.io/en/0.3.3/complete_example.html). The blast output will be a _PROKKA_XXXXXXX.out_ file which you can then give as an input to the script _ContigsCOGsDB.java_ in order to generate COG assignments for the protein sequences. Assuming that you are in the parent directory (CViewer), you can do this as follows:
```
$ java scripts/ContigsCOGsDB PROKKA_XXXXXXX.out
```
The above command will generate a TSV file of COGs named PROKKA_COGsDB.tsv placed in the _Output/_ folder, and should look like this:
```
$ head PROKKA_COGsDB.tsv
CDS_ID	COG	Name	Start	End
PROKKA_00008	COG0283	Cytidylate kinase [Nucleotide transport and metabolism]	1	189
PROKKA_00009	COG1302	Uncharacterized conserved protein YloU, alkaline shock protein (Asp23) family [Function unknown]	54	162
PROKKA_00010	COG1063	Threonine dehydrogenase or related Zn-dependent dehydrogenase [Amino acid transport and metabolism, General function prediction only]	3	339
PROKKA_00010	COG1064	D-arabinose 1-dehydrogenase, Zn-dependent alcohol dehydrogenase family [Carbohydrate transport and metabolism]	55	340
PROKKA_00011	COG1211	2-C-methyl-D-erythritol 4-phosphate cytidylyltransferase [Lipid transport and metabolism]	1	235
PROKKA_00011	COG0746	Molybdopterin-guanine dinucleotide biosynthesis protein A [Coenzyme transport and metabolism]	1	166
```

**Step 2: Importing the data into CViewer**

Once you have completed the above steps, then in a way similar to the one described in the previous section, click the **Open** button next to
* ```GenBank:``` to import the GenBank file (e.g. PROKKA.gbk)
* ```COGs:``` to import the file with all the identified COGs (e.g. PROKKA_COGsDB.tsv)

**Step 3: Extracting the annotation tracks**

Before displaying the annotation tracks, one has to extract them first from the GenBank file. This is conveniently automated in CViewer and one just needs to click and expand the ```Contig Functional Annotation``` section, check the ```Enable annotation``` _CheckBox_ and then press the ```Apply``` button. This will enable the software to extract the annotation tracks for all the data contigs from the GenBank file and save them into a condensed form for efficient processing afterwards. The data is saved in the _Output/_ folder under the name _AnnotatedContigsProkka.txt_. 

**Step 4: Visualising the annotation tracks**

Once we have successfully performed the above steps, one can visualise by interacting with the PCA plots and look for clustering of the contigs. When the user moves the cursor over the PCA plots, the position is captured and a label is then shown for the chosen contig with the name of the cluster that it belongs to. The user can then click on that contig to activate the accompanying panel where the data for the annotation tracks get displayed. When this is done, the CDS regions for the given contig are displayed in cyan colour, enzyme and COGs information is represented in yellow and red colours, respectively, while the labels for putative products for the CDS regions can also be replaced with the genes labels. 

**Step 5: Selecting and extracting annotation data**

The COGs identified in the CDS blocks are highlighted by right-clicking and choosing from a ```pop-up``` menu. If there is a need to explore some specific CDS regions further, then a set of them can be selected manually by the user and extracted as a ```fasta``` file with their sequence details (say a user wants to analyse how a particular homologous gene differs between different species through third party tools). Moreover, a legend describing summary statistics for the given contig is displayed on the bottom-left part of the panel including length, number of CDS regions, genes and enzyme identifiers. Finally, integration with statistical utilities for protein sequences allows for assessment of the structural characteristics of the CDS regions to highlight potential **```antimicrobial resistant genes```**. 

The above features are illustrated in the video below:
[![animation](https://user-images.githubusercontent.com/30604050/88488271-f8cbe980-cf83-11ea-87f5-02cbfba904d9.png)](https://www.youtube.com/watch?v=ZNXTTollJ7I)


## KEGG metabolic pathways
PROKKA software also assigns enzyme identifiers to the protein coding regions. The enzyme identifiers are represented by the [Enzyme Commission numbers](http://www.chem.qmw.ac.uk/iubmb/enzyme/) and describe a numerical classification scheme for the enzymes based on the chemical reactions they catalyse. In ```CViewer```, we have used these numbers to assign them to molecular level functions by utilizing the **[Kyoto encyclopaedia of genes and genomes (KEGG) Orthology database](https://www.genome.jp/kegg/ko.html)**. The database consists of a collection of functional orthologs or KEGG Ortholog (KO) groups, that each enzyme is part of, and correspond to the nodes (boxes) of the KEGG pathway maps. Using these KO groups and abundance data, the complete collection of molecular pathways for metabolism from KEGG database can be visualized in the software with annotation for the mappable enzymes within the pathway across all the samples. In addition to providing a visual cue for inspecting the pathways for a given sample, the software also allows for enrichment analysis of metabolic pathways between samples through downstream statistics, i.e., through Kruskal-Wallis test with appropriate adjustment for multiple comparisons (described in following section).

**Step 1: Extracting the KEGG Ortholog groups**

Before we start visualising the KEGG maps, we first need to extract the KO groups from the enzymes identified in the dataset. The matched KOs are first converted into presence/absence table for the contigs and combined next with the coverage profiles of that contigs to create microbial pathway abundances across the samples. This process is automated in CViewer and the only thing then required is to extract the data by pushing the Apply button right next to the label ```Extract KEGG Ortholog groups```. The software will then generate a ```CSV``` file describing all the KOs that were obtained for all the contigs across the samples and will place this file in the *Output/* folder.

**Step 2: Visualising the KEGG maps**

Once we have completed the above step, the complete collection of molecular pathways for metabolism from KEGG database can then be visualised and annotated in the software. This includes 11 key metabolic pathways with each of them describing multiple metabolic processes which can be chosen from a _drop-down_ list and visualised by pressing the ```Apply``` button located right next to the list. In each case, the metabolic maps are generated on a sample by sample basis i.e., one can visualise a map and inspect the KOs that are enriched in a sample of choice. To select a specific sample to inspect, one can use the same slider used for visualising the coverage profiles of the contigs through PCA plots that is provided in the ```Principal Component Analysis``` section. These steps are also demonstrated in the animation below:

[![animation](https://user-images.githubusercontent.com/30604050/88488454-151c5600-cf85-11ea-8865-c9f32f40e30b.png)](https://www.youtube.com/watch?v=21Ky1zD2mwU&feature=youtu.be)

## Phylogenetic diversity
Community phylogenetic alpha diversity analysis is also possible with ```CViewer``` software to explore environmental pressures and give an indication of stochasticity. We have considered two  metrics which, have been been suggested by Webb[[4]](#4) and look for phylogenetic clustering or overdispersion, followed by employing a randomisation procedure called statistical effect size (SES) to calculate the _Net Relatedness Index_ (NRI) and _Nearest Taxon Index_ (NTI). 

**Step 1: Extracting a gene from the metagenomics contigs**

Before extracting the gene data, we need to **ensure that the output from the CONCOCT (PCA file, coverage and clustering files; see Section: Visualising the metagenomics contigs) and PROKKA (GenBank file) pipelines are imported into the software**. Then, in order to perform phylogenetic analysis with CViewer, one first needs to generate a phylogenetic tree. To do this, the system provides a set of essential genes[[5]](#5) that can be extracted from the contigs for creating a phylogenetic tree with support for manual selection as well. After the gene has been specified, the gene sequences are extracted for each cluster of the metagenomics dataset and is saved in a ```.fasta``` file in the *Output/* folder, along with the average coverages of the clusters across all the samples. 

To illustrate the process described above, we will use the output from CONCOCT that was produced from the metagenomics profiling of 144 feacal extracts collected from Crohn's Disease (CD) patients undergoing dietary treatment with Exclusive Enteral Nutrition (EEN), Healthy Crohn's disease relatives (CDR), Healthy children (H) and Healthy adults (HC_A)[[6]](#6), [[7]](#7). Metagenomics reads were assembled into 636,062 contigs which, were subsequently grouped into 977 genomic bins through CONCOCT. The following video demonstrates how the gene data is extracted from the above dataset using CViewer. Please note that depending on the dataset size this process might require a few minutes to complete (approx. 12mins on a MacBook Pro with a 2.6 GHz Intel Core i5 processor):

[![animation](https://user-images.githubusercontent.com/30604050/89720701-1a3dc400-d9cd-11ea-87d0-fe84a284bd43.png)](https://www.youtube.com/watch?v=xD9wu2yXEOw)


**Step 2: Creating the phylogenetic tree**

Once we have extracted the sequences in a ```.fasta``` file, we can use the third party tools to generate the phylogenetic tree, such as [MUSCLE](http://www.drive5.com/muscle/) for performing the alignment and [FastTree](http://www.microbesonline.org/fasttree/) for producing the tree. Once this is accomplished, the obtained tree and the table describing the abundances across the dataset clusters (generated from the previous step) are uploaded to the software, to calculate the ***```NRI```*** and ***```NTI```*** indices. 

**Step 3: Phylogenetic alpha diversity analysis**

Using the tree and the abundance data obtained from the previous steps, community phylogenetic structure whether clustered or over-dispersed can be explored through the implemented phylogenetic α-diversity indices. The user can then choose between two different ```null models``` to generate the null communities. These include randomizations within samples while maintaining the species richness or within species while maintaining species frequency. In addition, it is possible to choose from a weighted (where the NRI and NTI metrics for each species are weighted by species abundance) or unweighted approach when performing the analysis. This is described below:

The data are first imported by clicking next to:
* ```Import tree:``` for the phylogenetic tree (e.g. alaS.tre)
* ```Import table:``` for the abundance data (e.g. alaS.csv)
* ```Import metadata:``` for the file with the associated metadata (e.g. CICRA_project_metadata.csv)

Then, the analysis can be configured based on the following parameters:
* ```Group by:``` the labels for the dataset samples as described in the metadata file (e.g. CICRA_project_metadata.csv)
* ```Maintain:``` whether randomizations should be done within samples while maintaining the species richness or within species while maintaining species frequency
* ```NRI```, ```NTI``` tickboxs: whether NRI or NTI should be estimated
* ```weighted``` tickbox: whether the NRI and NTI metrics for each species should be weighted by species abundance
* ```runs```: number of randomizations
* ```Show tree```: whether the phylogenetic tree should be displayed
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

The above steps are demonstrated in the following video:
[![animation](https://user-images.githubusercontent.com/30604050/89720708-24f85900-d9cd-11ea-8735-aef8cd65e3f4.png)](https://www.youtube.com/watch?v=ql0rcEjPeqs)


## Data analyses

CViewer, in addition to providing the methods for the visualisation and exploration of the features of contigs in the context of sample coverages, also supports a comprehensive set of multivariate statistical algorithms to allow exploratory as well as hypothesis driven analyses. The provided tools have been simplified and implemented in the form of intuitive workflows tailored to the need of non-expert users. For convenience, all these techniques have been categorised and included in the ```Methods``` tab of the software.

**Step 1: Preparing the data**

Almost every statistical procedure implemented in CViewer requires an abundance table and a metadata file providing the categorical labellings of the samples as well as extrinsic physico-chemical/clinical parameters. You can use ```CViewer``` to obtain an abundance table from the metagenomics contigs, by using the coverage information of the contigs across the dataset samples and the clustering file that was generated through the CONCOCT pipeline. To do this, the software provides a specific method that is found in the ```Main``` tab and is called ```Extract Table```. To generate the table, all you have to do is to press the ```Apply``` button right next to the label ```Extract average coverages```. This will generate a table describing the average coverages of the clusters which is subsequently and automatically put in the *Output/* folder of your current directory. You can also choose to label the dataset clusters with the taxonomic labels that were assigned to, if these are availabe, and which were imported into the software during the previous steps (see Section: Visualising the metagenomics contigs).


**Step 2: Importing the data into the software**

Once all the required information is uploaded, in principle, one can start the analysis without the need to import anything further. To import, you will need to navigate yourself to the ```Methods``` tab and then to the ```Input Data``` section of that tab. For demonstration purposes, in this tutorial, we have included an abundance table describing the average coverages across 977 clusters (obtained through CONCOCT software) and 144 samples collected from Crohn's Disease (CD) patients undergoing dietary treatment with Exclusive Enteral Nutrition, Healthy Crohn's disease relatives and Healthy individuals [[6]](#6), [[7]](#7). To import these data into the software, please do as follows:

In a way similar to the one described in previous sections, click the **Open** button next to
* ```Table:``` to import your abundance table (e.g. CICRA_average_coverages.csv)
* ```Metadata:``` to import the file with the associated metadata (e.g. CICRA_project_metadata.csv)

**Step 3: Normalization Techniques**

The tool offers a variety of popular normalization techniques to reduce systematic variation in data. These techniques include:

- *Relative Transformation*
- *Log Transformation* (Natural or Base 2 Logarithm)
- *Log-relative Transformation*
- *TSS+CLR Transformation*
- *Pareto Scaling* (Mainly for Metabolomics Data)
  
Users can select their desired normalization method from the dropdown list provided next to the data upload button. Each normalization technique serves a specific purpose and can be chosen based on the nature of the data and the analysis requirements.

**TSS+CLR Transformation**:

The TSS+CLR (Total Sum Scaling + Centered Log-Ratio) transformation is a specialized normalization technique designed for handling compositional data. Compositional data represent relative proportions of different parts within a whole and are commonly encountered in fields such as microbiome analysis.

**Total Sum Scaling (TSS)**:

TSS involves scaling each sample so that the total sum of all parts equals a constant value, typically 1. This ensures that the compositional nature of the data is preserved.
Centered Log-Ratio (CLR):

CLR transformation addresses the issue of compositional data being constrained to a simplex by transforming it into unconstrained space using logarithms. Additionally, CLR centers the data around the geometric mean, reducing the impact of outliers and improving the interpretability of the transformed data.

**Offset for Zeros**:

An offset of 1 is added by default to deal with zeros in the data, preventing mathematical errors during the transformation process.

## Alpha diversity
After we have performed the above steps, we can start exploring our data in CViewer. The tool allows alpha diversity analyses and considers a number of popular indices that can be used to assess the microbial diversities within a sample. We have considered the Shannon's (H'), Simpson’s diversity (D1) and its inverse (D2) which account for species richness and abundance. To measure how similar the distributions of species in a community are to each other, the tool provides the Pielou’s evenness. Finally, the relative proportions of the most dominant taxa for a given community dataset can also be explored. In each case, the generated plots can be **customised (e.g. change X-Y axis font size, range etc.) by right-clicking** on the plot and can be extracted into a ```PNG``` file for publication or further research. In the given menu you will find the options described below:

* ```Index:``` a drop-down menu with the available a-diversity methods (e.g. Shannon, Simpson etc.)
* ```Group by:``` the labels for the dataset samples as described in the metadata file (e.g. CICRA_project_metadata.csv)
* ```Taxa number:``` field activated when one wants to explore the ***Relative proportions*** of the most abundant clusters/species in a groups of samples and used for showing the number of the ***X*** top clusters/species in each group that one wishes to inspect.
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

See the animation below for more details:
[![animation](https://user-images.githubusercontent.com/30604050/89131606-fcbcb600-d505-11ea-89e5-dff538ceb576.png)](https://www.youtube.com/watch?v=qj6zGOsl_24)


## Beta diversity

### 1. Principal Component Analysis
Principal component analysis is one of the most popular among the existing dimensionality reduction techniques. The goal here is to find the best summary of the data using a small number of principal components (PCs). We have considered the Nonlinear Iterative Partial Least Squares (NIPALS)[[8]](#8) for computing PCA, a popular approach in multivariate data analysis that performs well with large datasets. One can explore up to 8 components of the PCA analysis using CViewer, where if one is interested just in the first few dimensions, the number of components to be calculated can then be specified manually to improve also time efficiency. The different components can be inspected by using the provided slider and the percentage of variance explained in each dimension is given in the resulting plot. In the menu you will find the following options:

* A ```CheckBox``` with the label ***All dimensions (max 8)***: default option used when the user has not specified the number of PCA components to be returned
* ```Num:``` field used to specify, if desired, the number of PCA components to be calculated
* ```Group by:``` the labels for the dataset samples as described in the metadata file (e.g. CICRA_project_metadata.csv)
* ```Components:``` slider-bar giving a means to shift between PCA components
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

In this example, we will use the same steps that were described in the previous section for Alpha diversity analysis for uploading the data and the same dataset. In this case, however, we will first normalise our data before generating the PCA plot. To do that, one can use the ```drop-down``` menu provided in the tool and described in section **Step 2: Normalising the abundance data**.

The video below provides an illustration of PCA analysis in CViewer:
[![animation](https://user-images.githubusercontent.com/30604050/89131611-01816a00-d506-11ea-88d7-5c807744a1d1.png)](https://www.youtube.com/watch?v=gznFDEfKz5M)

### 2. Multidimensional Scaling
Similar to PCA, Multidimensional scaling (MDS), also known as Principal Coordinate Analysis (PCoA), is another dimensionality reduction technique that attempts to represent the (dis)similarity between a set of objects in a reduced space, based on their pair-wise dissimilarities of samples given in the form of a distance matrix. The distance matrices can be computed in ```CViewer``` by choosing from some of the most commonly used measures, such as Bray-Curtis, Jaccard and Euclidean distance. The percentage of variance explained in each dimension is given in the resulting plot. The following options are provided for analysis using MDS:

* ```Distance:``` a drop-down menu with the available distances methods (Euclidean, Bray-Curtis and Jaccard)
* ```Group by:``` a drop-down menu with the columns containing labels for the dataset samples as described in the metadata file (e.g. CICRA_project_metadata.csv)
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

Following from the previous demonstration, in this example, we will assume that you have already imported the required data, i.e. abundance table and metadata, into the software. Then, the video belows shows how MDS is performed in CViewer:
[![animation](https://user-images.githubusercontent.com/30604050/89131612-02b29700-d506-11ea-93ad-e30fd18fbd9f.png)](https://www.youtube.com/watch?v=qj6zGOsl_24)

### Fuzzy Set Ordination
Fuzzy set ordination is used to test effects of pertubation in environmental variables to community structure (http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological/FSO_notes.pdf). For each of the specified variables, a fuzzy set ordination is calculated and the correlation between the original variable and the fuzzy set is reported. The significance of a particular variable is assessed by comparing a specified threshold p-value and the probability of obtaining a correlation between the data and fuzzy set.

In fuzzy set ordination, samples are assigned gradual membership values (fuzzy) ranging from 0 to 1. To achieve this, FSO does not use the raw community data, but rather a similarity (or distance) matrix, which is calculated prior to the statistic. The software supports three different similarity indices to calculate the similarity matrix, namely, the Baroni-Urbani & Buser, Horn and Yule indices. The results are visualised by producing a plot of fuzzy set against original values which is annotated with a correlation between them and a significance label. In the given menu, you will find the following options:

* ```Experimental variables:``` a ```CSV``` table describing the available experimental/environmental variables (Euclidean, Bray-Curtis and Jaccard) such as pH, Temperature, Calprotectin levels etc.
* ```Similarity index:``` a drop-down menu with the available similarity indices (Baroni-Urbani & Buser, Horn and Yule)
* ```Variable:``` field used to specify the environmetal variable of interest
* ```Group by:``` the labels for the dataset samples as described in the metadata file (e.g. CICRA_project_metadata.csv)
* ```Permutations:``` the number of permutations to be used for deriving the ***p-values*** (typically 1000)
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

In this example, we will assume again that you have already imported the required data, i.e. abundance table and metadata, into the software. Let's say that we are then interested in exploring if e.g. there is an association between the community composition of the gut microbiome of Crohn's disease patients with calprotectin levels or with the levels of SCFAs and how this differs between patients who achieved or not clinical remission at the end of EEN treatment (point D). Then, the video belows shows how FSO is performed for inspecting this in CViewer:
[![animation](https://user-images.githubusercontent.com/30604050/89131615-03e3c400-d506-11ea-8d08-b388771c1c09.png)](https://www.youtube.com/watch?v=_028HSzxUdw)

### Permutational Multivariate Analysis of Variance
Given a community dataset and a set of predictor physico-chemical variables, PERMANOVA can be used to provide information about the percentage of variation (R2) explained by the given predictors and the significance associated with it (P-value). The results of the PERMANOVA analysis are exported into a summary file in the *Output/* folder containing sources of variation, degrees of freedom, sequential sums of squares, mean squares, K statistics, partial R-squared and P-values, based on N permutations. In addition, the percentage (%) of variation explained (R2) by different predictor variables with annotation for significance can be visualised in the software in the form of a Bar Plot or a Pie chart. CViewer provides the following options for using PERMANOVA:

* A ```text-field```**```(Y~)```** used to provide the formula with the independent variables to be tested (e.g. Y ~ Group + pH + Calprotectin)
* ```Distance:``` a drop-down menu with the available distances methods (Euclidean, Bray-Curtis and Jaccard) that can be used to calculate pairwise distances 
* ```Style:``` results can be visualised in the form of a *Bar Plot* or a *Pie chart*
* ```Panel:``` which panelupper left/right and bottom left/right) should be used for plotting the results

This is demonstrated in the video below. We will assume again that you have already imported the required data, i.e. the abundance table and the metadata, into the software.:
[![animation](https://user-images.githubusercontent.com/30604050/89682919-6fa3a380-d8ef-11ea-9d62-f64af5ff7174.png)](https://www.youtube.com/watch?v=iga4uoq3U-w)


## Differential abundance

### Kruskal-Wallis H
Differential abundance analysis can be valuable when investigating for features (e.g. species/metabolic pathways) that discriminate between multiple treatment groups. ```CViewer``` utilizes the Kruskal-Wallis H statistic for this analysis. To enhance the credibility of the analysis amidst the challenge of multiple comparisons, CViewer incorporates the **Benjamini-Hochberg** correction method. This adjustment is crucial as it controls the **false discovery rate (FDR)**, which is the expected proportion of incorrect rejections among all rejected hypotheses. By adjusting P-values to control the proportion of type I errors (false positives) in the context of multiple comparisons, CViewer makes it less likely for users to mistakenly identify significant differences when there are none. This approach addresses the concern of false discovery in datasets with no real biological signal. In addition, pairwise significances between the groups can be explored with the Dunn’s post hoc procedure. The hypotheses to be tested can be uploaded in the software in the form of a .csv file containing hypotheses as headers along with additional metadata on the samples describing the treatment groups. One can then easily navigate through the list of available hypotheses and choose the one to be tested with the Kruskal-Wallis statistic. Finally, significant features are visualised in the form of box plots with appropriate colouring to indicate the treatment groups along with the P-values to report significance. The test provides the following options:

* ```P-Value:``` the P-value cut-off of the test
* ```Objects:``` how many objects or groups we would like to viusalise in the plot from the total number of the returned features that were found significant 
* ```Post-hoc test:``` whether a post-hoc test should be performed with the Kruskal-Wallis test (based on Dunn’s post hoc procedure)
* ```Adjust P-Values:``` whether P-Values should be adjusted to account for multiple comparisons. Enables the adjustment of P-Values using the **Benjamini-Hochberg** correction to account for the risk of false discoveries due to multiple comparisons.
* ```Group by:``` the labels for the dataset samples as described in the metadata file (e.g. CICRA_project_metadata.csv)
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

The example below illustrates how the above options are used in CViewer in order to perform differential abundance analysis. We will assume that you have already imported the required data, i.e. the abundance table and the metadata. Now, let's say that we are interested in exploring, e.g., whether any species (CONCOCT clusters) in the gut microbiome of CD patients changes significantly in abundance during the EEN treamtent, and how this differs from the healthy baseline:

[![animation](https://user-images.githubusercontent.com/30604050/89682931-75998480-d8ef-11ea-83c8-9548bcccf567.png)](https://www.youtube.com/watch?v=dTCAzbKhezg)

### Friedman Test
CViewer's analytical capabilities have been broadened with the implementation of the Friedman test, a non-parametric statistical test for detecting differences in treatments across multiple test attempts. This enhancement specifically addresses the complexities associated with repeated measures data, where the same subjects are exposed to different conditions over time or space. The test provides the following key features and implementation:

* ```Post-hoc test:``` To delve deeper into the specific treatment pairs that drive the observed differences, CViewer offers pairwise comparisons using Conover's F. This method adjusts for multiple comparisons and provides detailed insights into which treatments significantly differ from each other, enhancing the interpretability of the Friedman test results.

* ```Block by:``` In addition to the "Group by" functionality, which organizes samples based on treatment groups, the "Block by" option allows users to block the samples. This feature is essential for the Friedman test, enabling the analysis to account for the blocking factor that represents the subjects across which measurements are repeated. Blocking helps control for variability between subjects, ensuring that the differences detected are attributable to the treatments rather than inherent subject differences.

* ```Data Filtering and Information Message:``` CViewer performs internal data filtering to ensure the applicability of the Friedman test. This process includes checks for missing values and the appropriateness of the data structure for repeated measures analysis. Upon completion of this filtering, CViewer presents the user with an information message detailing the results of the filtering process.

By integrating the Friedman test alongside the Kruskal-Wallis test in CViewer, researchers can choose the appropriate statistical test based on the nature of their data—Kruskal-Wallis for independent samples and Friedman for repeated measures—ensuring analytical precision.


[![Picture 1](https://github.com/KociOrges/cviewer/assets/30604050/4910bbb3-16b6-4fad-aeb7-1e6aed496577)](https://youtu.be/IxbSKZ2NJVU?si=hPKpsTj7YZ1Z4HIm)


## Correlation
We have considered the Pearson’s product-moment coefficient to measure the degree of the association between two continuous variables, the Kendall’s tau coefficient to compare two quantities measured on at least an ordinal scale and the Spearman’s rank-order correlation as the nonparametric alternative of Pearson’s product-moment coefficient. Correlation analysis can be performed between two different data matrices (X and Y) that have the same set of samples in common. For example, after using CViewer to generate the abundance table with the average coverages of clusters across the samples of the input metagenomics dataset, one may want to explore how these clusters are associated with a number of environmental/experimental variables (e.g. pH, bacterial metabolites) that were measured on the same dataset samples and how these associations differentiate between the given treatments groups. 

Prior to the analysis, the input tables can be normalised, filtered and/or sorted according to a user-provided threshold. For a selected correlation coefficient, a full report can be displayed containing the correlation value (R) and the statistical significance (P-value and adjusted P-value) obtained for each pairwise comparison between the elements of the two data matrices. To do this, the software provides a number of options that can be used to configure the analysis and which are described below:

**Step 1: Input Data: Importing the data into the software**
* ```Table X:``` the X table, where rows describe samples and columns features/variables/species (e.g. CICRA_average_coverages.csv)
* ```Table Y:``` the Y table, where rows describe samples and columns features/variables/species (e.g. CICRA_bacterial_metabolites.csv)
* ```Normalise:``` normalisation can be performed individually for each matrix
* ```Metadata:``` the file with the associated metadata (e.g. CICRA_project_metadata.csv)

**Step 2: Options: Configuring the analysis**
* ```Method:``` the correlation coefficient of the test (Pearson’s product-moment, Kendall’s tau or Spearman’s rank-order correlation can be used)
* ```Adjust by:``` *P-Values* are adjusted using the *Benjamini-Hochbers's* procedure for multiple testing based on the ***X*** or ***Y*** table
* ```Filter X (>):``` keep the rows of ***table X*** for which the total sum is greater than that of the "filter X" value 
* ```Filter Y (>):``` keep the rows of ***table Y*** for which the total sum is greater than that of the "filter Y" value 
* ```Top (applied to X or Y table whether right next to Filter X or Filter Y):``` Return the first "Top" features from table ***X*** or ***Y*** depending on selection
* ```Sort (applied to X or Y table whether right next to Filter X or Filter Y):``` Sort the returned features of table ***X*** or ***Y*** (depending on selection) based on abundance in decreasing order
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

The video below shows how these options are used to perform Correlation analysis in CViewer. Let's assume that we interested in exploring for any underlying association between the metagenomics clusters or species of the WGS dataset and some additional variables that were measured on the same dataset samples:

[![animation](https://user-images.githubusercontent.com/30604050/89682930-73cfc100-d8ef-11ea-8083-3e5f202985fa.png)](https://www.youtube.com/watch?v=ohsiSI7r9lc)


## Integrated ‘omics analysis
Integrative analysis methods can provide better understanding and interpretation of the composite biological datasets and help elucidate the dynamics of the complex biological systems at a greater scale. ```CViewer``` implements three different techniques that have been proposed for simultaneous analysis of multiple omics datasets, namely Simultaneous Component Analysis with rotation to COmmon and DIstinctive components (DISCO-SCA)[[9]](#9), Joint and Individual Variation Explained (JIVE)[[10]](#10) and Two-way Orthogonal Partial Least Squares (O2PLS)[[11]](#11). All of these are helpful in providing a “global” view on the biological system of interest by decomposing the variability of the composite omics datasets into a joint variability or common structure. This represents the common themes underlying multiple omics datasets under study, and highlight individual variability or distinctive structure.

**Step 1: Collecting the data**

To illustrate how CViewer is useful for an integrated analysis of multi'omics datasets, we will use as previously, the metagenomics profile of Crohn's Disease patients undergoing dietary treatment with Exclusive Enteral Nutrition and healthy individuals[[6]](#6), [[7]](#7). In addition, as a second dataset, we will use a LC-MS metabolomics profile of 4,255 metabolites that were identified in 11 faecal extracts from eleven healthy children and to 54 faecal extracts from eleven children undergoing EEN for the treatment of active CD at timepoints before, during (15, 30, and 60 days), and after treatment[[12]](#12). To import these data into the software we will do as follows:

Click the **Open** button right next to:
* ```Table X:``` to import the first omics dataset ***X***, where rows describe samples and columns features/variables (e.g. CICRA_average_coverages.csv)
* ```Table Y:``` to import the second omics dataset ***Y***, where rows describe samples and columns features/variables (e.g. CICRA_metabolomics.csv)
* ```Normalise:``` normalisation can be performed individually for each dataset


**Step 2: Data pre-processing**

In an integrative or simultaneous analysis of multiple datasets, it is expected that the different blocks comprising the data are linked with a common set of entities. The tool supports an integrative analysis for two omics datasets that are linked by a common sample space, i.e. the measurements from the different omics data are obtained for the same set of samples (e.g. Crohn’s disease or Healthy individuals).

Prior to the analysis, it is useful that the data is pre-processed. As the data are generated from different omics technologies, they may be considerably different in size and scale. ```CViewer``` provides a number of pre-processing steps that can be applied to correct these differences. When the variables differ largely in scales (or abundance), one may consider ***centring*** the variables and/or ***scale*** them within each block to a sum of squares of one. In addition, ***weighting*** blocks together can be useful to avoid the effects of blocks having considerably different sizes.


**Step 3: Model selection**

Before we choose any of the given methods to perform component analysis, one must first provide the number of components that the composite dataset of omics data is expected to contain, along with their characterization as either common or distinctive. This information is required from each integrative approach (i.e., DISCO-SCA, JIVE and O2PLS) and is necessary so that the common and individual structures can be successfully identified and separated for each omics dataset. To estimate the optimal common and distinctive components, ```CViewer``` provides the Model Selection function. 

A component is characterised as common, if the ratios between the explained variances (SSQ) of each data block and its estimation (```SSQ_X/SSQ_(X_Y)``` and ```SSQ_Y/SSQ_(Y_X)```) is somewhere between 0.8 and 1.5. Subsequently, the process returns the optimal number of common components that were found to account for a sizeable amount of variance between the data blocks. Pre-processing of data should also be considered before performing the model selection analysis (see Step 1: Data pre-processing). For distinctive components, function allows the user to select the optimal number of distinctive components depending on the percentage of accumulated variance explained, the individual explained variance of each component, the absolute value of its variability or just a fixed number of components. This is done by configuring the following parameters:</p>

* ```PCA criteria:``` the PCA criteria for selection (%accum, single%, rel.abs, fixed.num)
* ```Group by:``` the labels for the dataset samples as described in the metadata file (e.g. CICRA_project_metadata.csv)
* ```R max:``` maximum common components
* ```Variance thres.:``` cumulative variance criteria for PCA selection for ```%accum``` or ```single%``` criteria
* ```Rel. variance thres.:``` relative variance criteria for ```fixed.num```
* ```PC components:``` fixed number of components for ```fixed.num```

See video below for more details:
[![animation](https://user-images.githubusercontent.com/30604050/89717262-852ad300-d9ac-11ea-9185-873a785bbd5e.png)](https://www.youtube.com/watch?v=c_HTw0MlYj4)

**Step 4: Configuring the analysis**

Integrated analysis of the given datasets can be performed by using either DISCO-SCA, JIVE or O2PLS. The individual data-blocks can be labelled and the number of common and distinctive components that were previously estimated using the Model Selection function can be provided as an input to the software. Pre-processing of the data is also supported and should be considered at this stage of the analysis (see Step 1: Data pre-processing). These features are provided through the following options:

* ```Method:``` the method of the integrative analysis (DISCO-SCA, JIVE or O2PLS)
* ```Names:``` the labels of the dasasets, one for each block of data, separated by comma (e.g. metagenome, metabolome)
* ```Common:``` the number of common components
* ```Distinctive:``` the number of distinctive components, one for each block of data, separated by comma


**Step 5: Plotting the results**

Once the above steps have been successfully performed, CViewer then provides a flexible framework for plotting the results. One can choose between common or individual results to be explored and can specify the 'omics dataset that one wants to inspect. The tool allows plotting of scores, loadings or both, for common and distinctive parts, while combined plots of both parts can also be produced. To do this, the **```Plots```** section provides the following options:

* ```What:``` either ```scores```, ```loadings``` or ```both```
* ```Type:``` either ```common```, ```individual``` or ```both```
* ```Block:``` which block to plot, either ```1``` or ```2```
* ```Combined:``` whether to make a plot of two components from one block, or components from different blocks
* ```Group by:``` the labels for the dataset samples as described in the metadata file (e.g. CICRA_project_metadata.csv)
* ```Comps.:``` the ***```x```***  and ***```y```*** components to plot for the chosen type and block. If ```Combined=true```, it indicates the component to plot for the first block and the component for the second block to plot together.
* ```Comps.:``` whether the labels for the dataset samples should be displayed on the plot
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

This is illustrated in the following video:
[![animation](https://user-images.githubusercontent.com/30604050/89717265-8956f080-d9ac-11ea-91cd-6f4650a75f57.png)](https://www.youtube.com/watch?v=H2HApGkVZnM)

Moreover, the loading contributions of each dataset to the common structure can be examined across the number of the available common components, where loading weights can be sorted according to importance based on absolute values. This is done by adjusting the following options:
* ```Contrib.:``` whether the loading contributions should be displayed
* ```Sort:``` whether the loading contributions should be sorted according to importance based on absolute values
* ```Comp.:``` the component of the common structure to plot
* ```Top:``` the number of the first ***```x```*** variables (e.g. species, metabolites) to plot
* ```Panel:``` which panel (upper left/right and bottom left/right) should be used for plotting the results

See animation below:
[![animation](https://user-images.githubusercontent.com/30604050/89717267-8b20b400-d9ac-11ea-9420-1351763eaf53.png)](https://www.youtube.com/watch?v=AL4s6EFxCY4)


## Charts
```CViewer``` provides also a separate tab for the visualisation of charts, with a Line Chart and a Box Plot Chart being supported by the software. This feature is primarily provided as an extra visualisation tool where it is possible to display information as a series of data points, particularly in the case of the Line Chart, and explore trends and patterns between pairs of variables (described from the X and Y axes). For each graph, it is possible to colour the data according to a specific grouping variable of interest. In addition, the Box Plot chart can represent every kind of interval observation using five different measures, i.e. the min, the max, the median and the lower and upper quartiles. This can be a very useful tool when ones need to display the distribution of data, identify outliers, and inspect if and how symmetrical or skewed the data are.

This is demonstrated in the following video:
[![animation](https://user-images.githubusercontent.com/30604050/89717609-194a6980-d9b0-11ea-8117-541ddf7d8523.png)](https://www.youtube.com/watch?v=LjxZyRe80xA)


## References
<a id="1">[1]</a>
Alneberg J, Bjarnason BS, De Bruijn I, et al. Binning metagenomic contigs by coverage and composition. Nat Methods. 2014;11(11):1144-1146. doi:10.1038/nmeth.3103
<br>
<br>
<a id="2">[2]</a>
Wood DE, Salzberg SL. Kraken: Ultrafast metagenomic sequence classification using exact alignments. Genome Biol. 2014;15(3). doi:10.1186/gb-2014-15-3-r46
<br>
<br>
<a id="3">[3]</a>
Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics. 2014;30(14):2068-2069. doi:10.1093/bioinformatics/btu153
<br>
<br>
<a id="4">[4]</a> 
Webb CO. Exploring the phylogenetic structure of ecological communities: An example for rain forest trees. Am Nat. 2000;156(2):145-155. doi:10.1086/303378
<br>
<br>
<a id="5">[5]</a> 
Campbell JH, O’Donoghue P, Campbell AG, et al. UGA is an additional glycine codon in uncultured SR1 bacteria from the human microbiota. Proc Natl Acad Sci U S A. 2013;110(14):5540-5545. doi:10.1073/pnas.1303090110
<br>
<br>
<a id="6">[6]</a> 
Quince C, Ijaz UZ, Loman N, et al. Extensive modulation of the fecal metagenome in children with Crohn’s disease during exclusive enteral nutrition. Am J Gastroenterol. 2015;110(12):1718-1729. doi:10.1038/ajg.2015.357
<br>
<br>
<a id="7">[7]</a> 
Gerasimidis K, Nikolaou CK, Edwards CA, McGrogan P. Serial fecal calprotectin changes in children with Crohn’s disease on treatment with exclusive enteral nutrition: associations with disease activity, treatment response, and prediction of a clinical relapse. J Clin Gastroenterol. 2011;45(3):234-239. doi:10.1097/MCG.0b013e3181f39af5
<br>
<br>
<a id="8">[8]</a> 
Wold H. Soft Modelling by Latent Variables: The Non-Linear Iterative Partial Least Squares (NIPALS) Approach. J Appl Probab. 1975;12(S1):117-142. doi:10.1017/s0021900200047604
<br>
<br>
<a id="9">[9]</a> 
Schouteden M, Van Deun K, Wilderjans TF, Van Mechelen I. DISCO-SCA. Behav Res Methods. 2014;46(2):576-587. doi:10.3758/s13428-013-0374-6
<br>
<br>
<a id="10">[10]</a> 
Lock EF, Hoadley KA, Marron JS, Nobel AB. Joint and individual variation explained (JIVE) for integrated analysis of multiple data types. Ann Appl Stat. 2013;7(1):523-542. doi:10.1214/12-AOAS597
<br>
<br>
<a id="11">[11]</a> 
Trygg J, Wold S. O2-PLS, a two-block (X-Y) latent variable regression (LVR) method with an integral OSC filter. In: Journal of Chemometrics. Vol 17. ; 2003:53- 64. doi:10.1002/cem.775
<br>
<br>
<a id="12">[12]</a> 
Alghamdi A, Gerasimidis K, Blackburn G, et al. Untargeted metabolomics of extracts from faecal samples demonstrates distinct differences between paediatric crohn’s disease patients and healthy controls but no significant changes resulting from exclusive enteral nutrition treatment. Metabolites. 2018;8(4). doi:10.3390/metabo8040082
