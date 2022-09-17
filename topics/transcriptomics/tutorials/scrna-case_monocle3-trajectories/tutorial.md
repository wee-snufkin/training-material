---
layout: tutorial_hands_on

title: 'Trajectory Analysis using Monocle3 '
subtopic: single-cell-CS
priority: 5
zenodo_link: 'https://zenodo.org/record/7078524'

questions:

objectives:

time_estimation: 1H

key_points:

requirements:
-
    type: "internal"
    topic_name: transcriptomics
    tutorials:
        - scrna-case_alevin
        - scrna-case_alevin-combine-datasets
        - scrna-case_basic-pipeline
tags:
- single-cell
- trajectory-analysis
- paper-replication

contributors:
- wee-snufkin
- nomadscientist

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

This tutorial is a follow-up to the 'Single-cell RNA-seq: Case Study' (find it [here](https://training.galaxyproject.org/training-material/topics/transcriptomics/)), we will use the same sample from the previous tutorials. If you haven’t done them yet, it’s highly recommended that you go through them to get an idea how to [prepare a single cell matrix](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-case_alevin/tutorial.html), [combine datasets](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-case_alevin-combine-datasets/tutorial.html) or [filter, plot and process scRNA-seq data](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-case_basic-pipeline/tutorial.html) to get the data in the form we’ll be working on today.

In this tutorial we will perform trajectory analysis using [monocle3](https://cole-trapnell-lab.github.io/monocle3/). You can find out more about the theory behind trajectory analysis here /slides/. We have already analysed the trajectory of our sample using ScanPy toolkit in another tutorial: [Trajectory Analysis using Python (Jupyter Notebook) in Galaxy](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-case_JUPYTER-trajectories/tutorial.html). However, trajectory analysis is quite sensitive and some methods work better for specific datasets. Now you can go through the same steps but using a different method to compare the results, usability and the final outcome! Sounds exciting, let’s dive into that! 

***tip: using tutorial mode***

## Get data
We still work on data from a mouse dataset of fetal growth restriction [Bacon et al. 2018](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-JUPYTER-trajectories/tutorial.html#Bacon2018) (see the study in Single Cell Expression Atlas [here](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and the project submission [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/)), making this case study even more comprehensive. 
Monocle3 works great with annotated data, so we will make use of our annotated AnnData object, generated in the previous [tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-seq-basic-pipeline/tutorial.html#neighborhood-graph). So you see - all the hard work of processing data was not in vain! We will also need a ‘clean’ expression matrix, extracted from the AnnData object just before we started the processing.
You can find both datasets in this [input history](https://humancellatlas.usegalaxy.eu/u/j.jakiela/h/monocle3-input-files) or download from Zenodo below.  


>### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the AnnData object from [Zenodo]({{ page.zenodo_link }})
>
>    ```
>    {{ page.zenodo_link }}/files/AnnData_before_processing.h5ad
>    {{ page.zenodo_link }}/files/Annotated_AnnData.h5ad
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Check that the datatype is `h5ad`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="h5ad" %}
>
{: .hands_on}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preparing the input files

## Extracting annotations

 As we want to use Monocle for the trajectory analysis, we will have to feed it with cell metadata, gene annotation and expression matrix files (in theory expression matrix alone could do, but then we wouldn’t have all those useful annotations that we were working on so hard!). In order to get those files, we will extract the gene and cell annotations from our AnnData object. 

 > ### {% icon question %} Questions
>
> How many lines do you expect to be in gene annotations and cell metadata files?
>
> > ### {% icon solution %} Solution
> >
> > If you click on the step with uploaded annotated AnnData file, you will see on a small preview that this object has 8605 observations and 15395 variables, so we expect to get a cell metadata file with 8605 lines and gene annotations file with 15395 lines (without headers of course!).
> >
> {: .solution}
>
{: .question}

> ### {% icon hands_on %} Hands-on: Extracting annotations
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Annotated_AnnData` 
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
> 2. Rename {% icon galaxy-pencil %} the observations annotation `Extracted cell annotations (obs)`
>
> 3. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Annotated_AnnData` 
>    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
>
> 4. Rename {% icon galaxy-pencil %} the annotation of variables `Extracted gene annotations (var)`
>
>
{: .hands_on}

Quick and easy, isn’t it? However, there are some minor changes that we have to make to our files first. 

## Cell metadata
Our current dataset is not just T-cells: as you might remember from the last tutorial, we identified a cluster of macrophages as well. This might be a problem, because the trajectory algorithm will try to find relationships between all the cells (even if they are not necessarily related!), and not only the T-cells that we are interested in. We need to remove those unwanted cell types to make the analysis more biologically relevant.

Manipulate AnnData tool allows you to filter observations or variables and that would be the easiest way to do it! However, we have to think ahead. As we want to use Monocle later on, we will have to provide it with cell metadata, gene annotation and expression matrix files anyway. This is why we extract the annotations first, make changes to them and finally we adjust the expression matrix to the filtered annotations. In that way, we’ll end up with three separate files, ready to be passed onto Monocle3.

 > ### {% icon question %} Questions
>
> Where is the information about cell types stored?
>
> > ### {% icon solution %} Solution
> >
> > We have already extracted the cell annotations file - in one of the columns you can find the information about cell type, assigned to each cell. 
> >
> {: .solution}
>
{: .question}

Click on `Extracted cell annotations (obs)` file to see a small preview window. This shows you that the column containing the cell types has number 22.  We’ll need that to filter out unwanted cell types!

> ### {% icon warning %} Check the column number!
> If you are working on a different dataset, the number of the ‘cell_type’ column might be different, so make sure you check it on a preview and use the correct number! 
{: .warning}

> ### {% icon hands_on %} Hands-on: Filter out macrophages
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `Extracted cell annotations (obs)`
>    - *"With following condition"*: `c22!='Macrophages'`
>    - *"Number of header lines to skip"*: `1`
>    - That’s it - our cell annotation file is ready for Monocle! Let’s rename it accordingly. 
> 2. **Rename** {% icon galaxy-pencil %} the output: `Cells input data for Monocle3`
>
>    > ### {% icon details %} Details: Parameters
>    >
>    > - `c22` means column no. 22 - that's the column with cell types, and it will be filtered for the macrophages
>    > - `!=` means 'not equal to' - we want to keep the cell types which ARE NOT macrophages
>    {: .details}
>
>    > ### {% icon tip %} Other unwanted cell types
>    >
>    > It might happen that during clustering you’ll find another cell type that you want to get rid of for the trajectory analysis. Then simply re-run this tool on already filtered file and change ‘Macrophages’ to another unwanted cell type.
>    {: .tip}
{: .hands_on}

## Gene annotations
Sometimes certain functionalities require specific indication where the data should be taken from. In case of Monocle3, to allow further genes analysis using one of its functions, it is essential that the names of the genes are stored in a column called ‘gene_short_name’. Therefore, we need to check what is the name of that column in our dataset. 

> ### {% icon question %} Questions
>
> 1. Where can you check the header of a column containing genes names?
> 2. What is the name of this column?
>
> > ### {% icon solution %} Solution
> >
> > 1. Our extracted gene annotations file! Either by clicking on the eye icon {% icon solution %} or having a look at the small preview window. 
> > 2. In our dataset the gene names are stored in a column called ‘Symbol’ - we need to change that!
> >
> {: .solution}
>
{: .question}

Let’s click on the `Extracted gene annotations (var)` file to see a small preview. We can see that the gene names are in the third column with a header ‘Symbol’. Keep that in mind - we’ll use that in a second!

> ### {% icon hands_on %} Hands-on: Changing the colname
>
> 1. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.2) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `Extracted gene annotations (var)` 
>    - *"using column"*: `c3`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Symbol`
>            - *"Replacement"*: `gene_short_name`
>     - Voila! That’s the gene input for Monocle! Just a quick rename...
> 2. **Rename** {% icon galaxy-pencil %} the output: `Genes input data for Monocle3`
>
{: .hands_on}

## Expression matrix
Last, but not least! I would even dare to say that last, but the most important! Actually the expression matrix contains all the values representing expression level of a particular gene in a cell. This is why in theory the expression matrix is the only input file required by Monocle3 - without annotation files the CDS data can still be generated - in fact it will be quite bare, but at least it could be processed. 

So, the values in the expression matrix are just some numbers. But do you remember that we have already done some processing such as normalisation and calculation of principal components on AnnData in the previous tutorial? That affected our expression matrix. Preprocessing is one of the steps in the Monocle3 workflow, so we want to make sure that the calculations are done on a ‘clean’ expression matrix. If we apply too many operations on our raw data, it will be too ‘deformed’ to be reliable. The point of the analysis is to use algorithms that make the enormous amount of data understandable in order to draw meaningful conclusions in accordance with biology. 

So how do we do that?
> ### {% icon question %} Questions
>
> 1. How many cells and genes are there in the `Anndata_before_processing` file? 
> 2. How many lines are there in `Cells input data for Monocle3`?
> 3. How many lines are there in `Genes input data for Monocle3`?
>
> > ### {% icon solution %} Solution
> > You can answer all the questions just by clicking on the given file and looking at the preview window.
> > 1. [n_obs x n_vars] = 31178 x 35734, so there are 31178 cells and 35734 genes.
> > 2. 8570 lines, including a header, which makes 8569 cells. 
> > 3. 15396 lines, including a header, which makes 15395 genes.
> >
> {: .solution}
>
{: .question}

As you can see, there are way more genes and cells in the unprocessed AnnData file, so the expression matrix is much bigger than we need it to be. If the genes and cells we prepared for Monocle3 are not the same as in the expression matrix, Monocle3 will crash. Therefore, we have to filter that big, clean matrix and adjust it to our already prepared genes and cells files. But first, let’s extract this matrix from the unprocessed AnnData object. 

> ### {% icon hands_on %} Hands-on: Extracting matrix
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `AnnData_before_processing` 
>    - *"What to inspect?"*: `The full data matrix`
> 2. **Rename** {% icon galaxy-pencil %} the output: `Unprocessed expression matrix`
>
{: .hands_on}

If you have a look at the preview of `Unprocessed expression matrix`, you’ll see that the first column contains the cell barcodes, while the first row - the gene IDs. We would like to keep only the values corresponding to the cells and genes that are included in `Cells input data for Monocle3` and `Genes input data for Monocle3`. How do we do it? First, we compare the cell barcodes from `Cells input data for Monocle3` to those in `Unprocessed expression matrix` and ask Galaxy to keep the values of the matrix for which the barcodes in both files are the same. Then, we’ll do the same for gene IDs. So we have to cut the first columns from `Cells input data for Monocle3` and `Genes input data for Monocle3` to be able to compare those columns side by side with the matrix file.

> ### {% icon hands_on %} Hands-on: Cutting out the columns
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `Cells input data for Monocle3`
> 2. **Rename** {% icon galaxy-pencil %} the output: `Cells IDs`
> 3. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `Genes input data for Monocle3`
> 4. **Rename** {% icon galaxy-pencil %} the output: `Genes IDs`
>
{: .hands_on}

> ### {% icon hands_on %} Hands-on:  Filter matrix values by cell barcodes
>
> 1. {% tool [Join two Datasets](join1) %} with the following parameters:
>    - {% icon param-file %} *"Join"*: `Cells IDs` 
>    - *"using column"*: `c1`
>    - {% icon param-file %} *"with"*: `Unprocessed expression matrix`
>    - *"and column"*: `c1`
>    - *"Keep lines of first input that do not join with second input"*: `Yes`
>    - *"Keep lines of first input that are incomplete"*: `Yes`
>    - *"Fill empty columns"*: `No`
>    - *"Keep the header lines"*: `Yes`
> 2. **Rename** {% icon galaxy-pencil %} the output: `Pre-filtered matrix (by cells)`
>
{: .hands_on}

Look at the preview of the output file. First of all, you can see that there are 8570 lines (8569 cells) instead of 31178 cells that were present in the matrix. That’s exactly what we wanted to achieve - now we have information for the T-cells that we had filtered. However, the step that we have already performed left us with the matrix whose first and second columns are the same - let’s get rid of one of those! 

> ### {% icon hands_on %} Hands-on: Remove duplicate column (cells IDs)
>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `Pre-filtered matrix (by cells)`
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `c1`
> 2. **Rename** {% icon galaxy-pencil %} the output: `Filtered matrix (by cells)`
>
{: .hands_on}

Now we will perform the same steps, but for gene IDs. But gene IDs are currently in the first row, so we need to transpose the matrix, and from there we can repeat the same steps as above, but for gene IDs of course. 

> ### {% icon hands_on %} Hands-on: Filter matrix by gene IDs
>
> 1. {% tool [Transpose](toolshed.g2.bx.psu.edu/repos/iuc/datamash_transpose/datamash_transpose/1.1.0+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `Filtered matrix (by cells)` 
>    - The matrix is now ready to be filtered by gene IDs!
> 2. {% tool [Join two Datasets](join1) %} with the following parameters:
>    - {% icon param-file %} *"Join"*: `Genes IDs`
>    - *"using column"*: `c1`
>    - {% icon param-file %} *"with"*: output of **Transpose** {% icon tool %}
>    - *"and column"*: `c1`
>    - *"Keep lines of first input that do not join with second input"*: `Yes`
>    - *"Keep lines of first input that are incomplete"*: `Yes`
>    - *"Fill empty columns"*: `No`
>    - *"Keep the header lines"*: `Yes`
> 3. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: output of **Join two Datasets** {% icon tool %}
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `c1`
>    -  Monocle3 requires that in the matrix rows are genes, and columns are cells - that is what we've got, so there is no need to transpose matrix again. The expression matrix is ready! Let's just rename it...
> 4. **Rename** {% icon galaxy-pencil %} the output: `Expression matrix for Monocle3`
>
{: .hands_on}

{% icon trophy %} Finally! We have prepared all the files to pass them onto the Monocle3 workflow!


# Monocle3 workflow

What will happen with those files that we have been preparing so far? Well, Monocle3 turns the expression matrix, cell and gene annotations into an object called cell_data_set (CDS), which holds single-cell expression data. 

> ### {% icon details %} Details: Input files
> 
> That’s what [Monocle3 documentation](https://cole-trapnell-lab.github.io/monocle3/docs/starting/) says about the required three input files:
>    - expression_matrix, a numeric matrix of expression values, where rows are genes, and columns are cells. Must have the same number of columns as the cell_metadata has rows and the same number of rows as the gene_metadata has rows.
>    - cell_metadata, a data frame, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
>    - gene_metadata, a data frame, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc. One of its columns should be named "gene_short_name", which represents the gene symbol or simple name (generally used for plotting) for each gene.
>
{: .details}

Here’s how the Monocle3 workflow looks like:

![Alternative text](../../images/image_name "Workflow provided by Monocle3 documentation")

We will follow those steps and see how it all works in practice. 


> ### {% icon hands_on %} Hands-on: Create CDS object
>
>    > ### {% icon details %} Details: Data format
>    >
>    > You can provide expression matrix as TSV, CSV, MTX or RDS file, while genes and cells metadata as TSV, CSV or RDS files. In our case all three files are tabular, so we will set the format to TSV.
>    {: .details}
> 1. {% tool [Monocle3 create](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_create/monocle3_create/0.1.4+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Expression matrix, genes as rows, cells as columns. Required input. Provide as TSV, CSV or RDS."*: `Expression matrix for Monocle3`
>    - *"Format of expression matrix"*: `TSV`
>    - {% icon param-file %} *"Per-cell annotation, optional. Row names must match the column names of the expression matrix. Provide as TSV, CSV or RDS."*: `Cells input data for Monocle3`
>    - *"Format of cell metadata"*: `TSV`
>    - {% icon param-file %} *"Per-gene annotation, optional. Row names must match the row names of the expression matrix. Provide as TSV, CSV or RDS."*: `Genes input data for Monocle3`
>    - *"Format of gene annotation"*: `TSV`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> What are the dimensions of the created CDS object?
>
> > ### {% icon solution %} Solution
> >
> > Just click on the performed step - on the preview you’ll see that the dimensions are 15395 x 8569 - so exactly as we predicted genes x cells! 
> >
> {: .solution}
{: .question}

## Pre-processing and dimensionality reduction

In Galaxy, there are currently 2 methods of initial dimensionality reduction which is included in the pre-processing step: principal component analysis (PCA) and latent semantic indexing (LSI). 
However, PCA is more commonly used, and it also allows us to perform further steps on CDS object, so we’ll use this method. There is one parameter here that has a great impact on how our analysis will look like, namely - the dimensionality of the initially reduced space. After many trials and errors, we were finally able to find the value that yielded the best results. You can have a look at the image below to see how different values affect the outcome - I can tell you now that we’ll go ahead with the value of **250**. Don’t worry, after a few more steps you’ll understand what all those colors mean and how to generate those plots. 

![Alternative text](../../images/image_name "Different outputs depending on the number of dimensions that the space is reduced to.")

> ### {% icon details %} Details: PCA & LSI
> 
> EXPLAIN PCA & LSI
>
{: .details}

> ### {% icon hands_on %} Hands-on: Pre-processing
>
> 1. {% tool [Monocle3 preprocess](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_preprocess/monocle3_preprocess/0.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 create** {% icon tool %}
>    - *"The dimensionality of the reduced space."*: `250`
 >
{: .hands_on}

Now it’s time for the proper dimensionality reduction so that instead of the initial thousands of dimensions, we can have only 2 and hence plot all the cells on one 2D graph. Again, there are several algorithms to do that: UMAP, tSNE, PCA and LSI (only possible when preprocess_method is set to 'LSI' as well), but due to the same reasons as above, we’ll use UMAP (most common + allows further operations). But I’ll let you see how the output from other algorithms look to convince you that **UMAP** is indeed the best in this case. 

![Alternative text](../../images/image_name "Different outputs depending on the algorithm of dimentionality reduction, applied to the output of the previous step (except LSI method which was called on LSI-preprocessed data).")

> ### {% icon hands_on %} Hands-on: Dimensionality reduction
>
> 1. {% tool [Monocle3 reduceDim](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_reducedim/monocle3_reduceDim/0.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: `output_rds` (output of **Monocle3 preprocess** {% icon tool %})
>
{: .hands_on}

## Plotting
 
Alright, now let's have a look at our output! Above you got a sneak peek of how the plot would look like, but now you’ll generate the plot on your own! 
Thanks to the fact that we provided Monocle3 with annotated data, we can now color the cells by any attribute that was in the cell metadata file! So, similarly to the previous tutorial, we’ll color them by cell type, genotype, batch and sex. At least for now. 


> ### {% icon hands_on %} Hands-on: Plotting
>
> 1. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 reduceDim** {% icon tool %}
>    - *"The cell attribute (e.g. the column of pData(cds)) to map to each cell's color, or one of {cluster, partition, pseudotime}."*: `cell_type`
> 2. Rename {% icon galaxy-pencil %} the output: `Cell type plot`
>
> 3. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 reduceDim** {% icon tool %}
>    - *"The cell attribute (e.g. the column of pData(cds)) to map to each cell's color, or one of {cluster, partition, pseudotime}."*: `genotype`
>    - *"If set, display the cell group names directly on the plot. Otherwise include a color legend on the side of the plot."*: {% icon history-share %} `No`
> 4. Rename {% icon galaxy-pencil %} the output: `Genotype plot`
>
> 5. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 reduceDim** {% icon tool %}
>    - *"The cell attribute (e.g. the column of pData(cds)) to map to each cell's color, or one of {cluster, partition, pseudotime}."*: `batch`
>    - *"If set, display the cell group names directly on the plot. Otherwise include a color legend on the side of the plot."*: {% icon history-share %} `No`
> 6. Rename {% icon galaxy-pencil %} the output: `Batch plot`
>
> 7. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 reduceDim** {% icon tool %}
>    - *"The cell attribute (e.g. the column of pData(cds)) to map to each cell's color, or one of {cluster, partition, pseudotime}."*: `sex`
>    - *"If set, display the cell group names directly on the plot. Otherwise include a color legend on the side of the plot."*: {% icon history-share %} `No`
> 8. Rename {% icon galaxy-pencil %} the output: `Sex plot`
>
{: .hands_on}


***comments, images***

## Clustering

Before inferring the trajectory, we have to group cells into clusters, which is an important step in identifying the cell types represented in your data. Monocle uses a technique called [community detection](https://doi.org/10.1038/s41598-019-41695-z)  to group cells. This approach was introduced by [Levine et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/) as part of the phenoGraph algorithm. Monocle also divides the cells into larger, more well separated groups called partitions, using a statistical test from [Alex Wolf et al](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1663-x), introduced as part of their [PAGA](https://github.com/theislab/paga) algorithm.

> ### {% icon details %} Details: Clusters vs partitions
> 
> Clusters are particularly useful while trying to assign cells to a certain type, because they are based on the similarity in gene expression. While inferring the trajectory, we will be analysing the relationships between clusters.
>Partitions are larger groups of cells, usually containing several clusters. Trajectory inference is performed only within one partition, so it is essential that all the cells that we want to analyse in pseudotime belong to the same partition. 
>
{: .details}


> ### {% icon hands_on %} Hands-on: Clustering 
>
> 1. {% tool [Monocle3 partition](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_partition/monocle3_partition/0.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 reduceDim** {% icon tool %}
>    - *"The q-value threshold used to determine the partition of cells."*: `1.0`
>    - The clusters and partitions are now stored in your CDS file. To see them, just plot the output, coloring the cells with the corresponding attributes.
>
> 2. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 partition** {% icon tool %}
>    - *"The cell attribute (e.g. the column of pData(cds)) to map to each cell's color, or one of {cluster, partition, pseudotime}."*: `partition`
> 3. Rename {% icon galaxy-pencil %} the output: `Partition plot`
>
> 4. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 partition** {% icon tool %}
>    - *"The cell attribute (e.g. the column of pData(cds)) to map to each cell's color, or one of {cluster, partition, pseudotime}."*: `cluster`
> 5. Rename {% icon galaxy-pencil %} the output: `Cluster plot`
>
>    > ### {% icon tip %} If the granularity of clusters is not satisfying...
>    >
>    > If you are not satisfied with the results of the standard igraph louvain clustering algorithm, you may set the `resolution` of clustering, which specifies the granularity of clusters. 
>    > ![Alternative text](../../images/image_name "On the left - clusters formed using standard igraph louvain clustering algorithm and on the right - clusters formed when the resolution was set to")
>    {: .tip}
>    > ### {% icon tip %} If the partition does not contain all cells of interest...
>    >
>    > Sometimes it might happen that cells are grouped into several partitions, while you want them all to be in just one in order to perform trajectory analysis on all of them. Then, you can try to increase the `q-value` threshold that is used to determine the partition of cells. 
>    > ![Alternative text](../../images/image_name "Describe the changes")
>    {: .tip}
>
{: .hands_on}

 
**compare clusters and cell types** - jpg? 
![Alternative text](../../images/image_name "Describe the changes")
 

## Gene expression
>We haven't looked at gene expression yet! This step is particularly important when working with data which is not annotated. Then, based on the expression of marker genes, you are able to identify which clusters correspond to which cell types. This is indeed what we did in the previous tutorial using scanpy. We can do the same using Monocle3! Since we work on annotated data, we can directly check if the expressed genes actually correspond to the previously assigned cell types. If they do, that’s great - if two different methods are consistent, that gives us more confidence that our results are valid. 
>Below is the table that we used in the previous tutorial to identify the cell types.

| Marker | Cell type |
|--------------------|
| Il2ra    | Double negative (early T-cell)    |
| Cd8b1, Cd8a, Cd4    | Double positive (middle T-cell)|
| Cd8b1, Cd8a, Cd4 - high | Double positive (late middle T-cell)|
| Itm2a    | Mature T-cell |
| Aif1    | Macrophages    |
| Hba-a1    | RBC    |

> ### {% icon hands_on %} Hands-on: Gene expression
>
> 1. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 partition** {% icon tool %}
>    - *"The cell attribute (e.g. the column of pData(cds)) to map to each cell's color, or one of {cluster, partition, pseudotime}."*: `cell_type`
>    - *"A list of gene IDs/short names to plot."*: `Il2ra,Cd8b1,Cd8a,Cd4,Itm2a,Aif1,Hba-a1`
>
{: .hands_on}

Let’s look at the expression of those genes using Monocle3 and compare them with the `Cell type plot`
![Alternative text](../../images/image_name "Describe the changes")

**analysis**

Here we used a priori knowledge regarding the marker genes. If we wanted to approach this problem in an unsupervised manner, we could use Monocle to tell us what would be the top marker genes in each group of cells. This is very useful if we don’t know the type of the cells in a specific cluster and we want to identify it, based on common marker genes. Or when we want to find other marker genes than those currently known.

> ### {% icon hands_on %} Hands-on: Top marker genes
>
> 1. {% tool [Monocle3 top markers](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_topmarkers/monocle3_topmarkers/0.1.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Object"*: output of **Monocle3 partition** {% icon tool %}
>    - *"Group cell by"*: `cell_type`
> 2. Rename {% icon galaxy-pencil %} the tabular output: `Top markers table`
> 3. Rename {% icon galaxy-pencil %} the pdf output: `Top markers plot`
>
{: .hands_on}

**analysis**

> ### {% icon question %} Questions
>
> If I cluster cells that are not annotated, can I assign clusters to cell type based on gene expression using Monocle3?
>
> > ### {% icon solution %} Solution
> >
> > Of course you can! That’s the point of clustering and gene expression analysis. However currently this function hasn’t been turned into a Galaxy tool yet, so in order to do so, you have to use the piece of code in R which performs this annotation. 
> >
> {: .solution}
>
{: .question}

There is also one more tool that allows you to get a better insight into gene expression, namely - identifying differentially expressed genes along the inferred trajectory. But in order to do that, we have to infer that trajectory first!

## Learn the trajectory graph

We’re getting closer and closer! The next step is to learn the trajectory graph, which means to fit a principal graph within each partition. In that way, we’ll ‘connect’ the existing clusters by creating a path between them.

> ### {% icon hands_on %} Hands-on: Learn graph
>
> 1. {% tool [Monocle3 learnGraph](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_learngraph/monocle3_learnGraph/0.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 partition** {% icon tool %}
>    - Again, the graph is now stored in your CDS file. To see it, just plot the output, you can color the cells by any attribute that you want. We'll use cell types to see how they are connected.
> 2. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 learnGraph** {% icon tool %}
>    - *"The cell attribute (e.g. the column of pData(cds)) to map to each cell's color, or one of {cluster, partition, pseudotime}."*: `cell_type`
>
{: .hands_on}

**analysis**

## Pseudotime analysis

Finally it's time to see our cells in pseudotime! We already learned trajectory, now we only have to order cells along it. Monocle3 requires information where to start ordering the cells, so we need to provide it with this information. We annotated early T-cells as double negative (DN), so those will be our root cells! 

> ### {% icon details %} Details: Pseudotime
> 
> To infer trajectories, we need data from cells at different points along a path of differentiation. This inferred temporal dimension is known as pseudotime. Pseudotime measures the cells’ progress through the transition. 
[read more](https://doi.org/10.1093%2Fbioinformatics%2Fbtw372)
>
{: .details}

> ### {% icon hands_on %} Hands-on: Ordering the cells along trajectory
>
> 1. {% tool [Monocle3 orderCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_ordercells/monocle3_orderCells/0.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 learnGraph** {% icon tool %}
>    - *"The cell phenotype (column in pdata) used to identify root principal nodes."*: `cell_type`
>    - *"The value in the cell phenotype column used to extract root nodes."*: `DN`
>    - Alright - we were waiting for this plot the whole tutorial: once we have the cells ordered, we can finally color them by pseudotime!
>
> 2. {% tool [Monocle3 plotCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_plotcells/monocle3_plotCells/0.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 orderCells** {% icon tool %}
> 3. Rename {% icon galaxy-pencil %} the output: `Pseudotime plot`
>
{: .hands_on}

**analysis**

Once the trajectory has been inferred, you might want to return to the gene expression analysis and dive into that in more depth. Here is a powerful tool that would give you even more information about the genes. 


> ### {% icon hands_on %} Hands-on: Differentially expressed genes
>
> 1. {% tool [Monocle3 diffExp](toolshed.g2.bx.psu.edu/repos/ebi-gxa/monocle3_diffexp/monocle3_diffExp/0.1.4+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in RDS format"*: output of **Monocle3 orderCells** {% icon tool %}
> 2. Rename {% icon galaxy-pencil %} the output: `Differential gene expression table`
>
{: .hands_on}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
