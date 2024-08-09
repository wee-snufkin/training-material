---
layout: tutorial_hands_on
title: Predicting EI+ mass spectra with QCxMS
zenodo_link: ' https://zenodo.org/record/13259853 '
level: Intermediate

questions:
- Can I predict QC-based mass spectra starting from SMILES?
- Do I need HPC environment to run the predictions? 
- Can I take into account different conformers?

objectives: 
- To prepare SDF files for downstream analysis, starting from SMILES. 
- To generate conformers and optimise them using semi-empirical methods.
- To produce simulated mass spectra for a given molecule in MSP format.  

time_estimation: 1H

key_points:
- Galaxy provides access to HPC resources and hence allows the use of semi-empirical methods and in silico prediction of QC-based mass spectra.
- 

tags:

contributions:
  authorship:
    - wee-snufkin
    - hechth
    - xtrojak
    - maximskorik

requirements :
  - type: "internal"
    topic_name: metabolomics
    tutorials: 
      - lcms-preprocessing

follow_up_training:
  -
    type: "internal"
    topic_name: 
---



> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


We will start with a table with the first column being molecule names and the second one – corresponding SMILES. 

INFO ABOUT SMILES

## Get data and pre-process
You have three options for uploading the data. The first two - importing via history and Zenodo link will give a file specific to this tutorial, while the last one – “Paste data” uploader gives you more flexibility in terms of the compounds you would like to test with this workflow. 

> <hands-on-title>Option 1: Data upload - Import history</hands-on-title>
>
> 1. You can simply import [history](https://usegalaxy.eu/u/j.jakiela/h/input-file-end-to-end-ei-mass-spectra-prediction-workflow-using-qcxms) with the input table. 
>
>    {% snippet faqs/galaxy/histories_import.md %}
>
> 2. **Rename** {% icon galaxy-pencil %} the history to your name of choice.
>
{: .hands_on}

><hands-on-title>Option 2: Data upload - Add to history via Zenodo</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the input table from [Zenodo]({{ page.zenodo_link }})
>
>    ```
>    https://zenodo.org/records/13259853/files/qcxms_prediction_input.tabular
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

> <hands-on-title> Option 3: Data Upload  - paste data </hands-on-title>
> 
> 1. Create a new history for this tutorial
> 2. * Click {% icon galaxy-upload %} **Upload Data** at the top of the tool panel
> * Select {% icon galaxy-wf-edit %} **Paste/Fetch Data** at the bottom
> * Paste the contents into the text field, separated by space. First, enter the name of the molecule, then its SMILES. Please note that we are not using headers here. For this tutorial, we’ll use the example of ethanol and ethylene, but feel free to use your own examples. 
> ```
> ethanol CCO
> ethylene C=C
> ```
> * Change **Type** from "Auto-detect" to `tabular`
> * Find the gear symbol ({% icon galaxy-gear %}), and select only **Convert spaces to tabs**
> * Press **Start** and **Close** the window
>
> 3. You can then rename the dataset as you wish.
> 4. Check that the datatype is *tabular*.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

Once your dataset has been uploaded, we can do some simple pre-processing to prepare the file for downstream analysis. Let’s start with ‘cutting’ the table into two columns – one with SMILES, the other with the molecule name. 

> <hands-on-title> Cutting out name column </hands-on-title>
>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/9.3+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: the tabular file in your history with the compound name and SMILES
>    - *"Operation"*: `Keep`
>    - *"Cut by"*: `fields`
>    - *"Delimited by"*: `Tab`
>    - *"Is there a header for the data’s columns"*: `No`
>    - *"List of Fields"*: `1`
>
{: .hands_on}

> <hands-on-title>  </hands-on-title>
>
> 1. {% tool [Split file to dataset collection](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_to_collection/split_file_to_collection/0.5.2) %} with the following parameters:
>    - {% icon param-file %} *"Select the file type to split"*: `Tabular`
>    - *"Tabular file to split"*: ``
>    - *"Number of header lines to transfer to new files"*: `0`
>    - *"Split by row or by a column?"*: `By row`
>    - *"Specify number of output files or number of records per file?"*: `Number of records per file (‘chunk mode’)`
>    - *"Chunk size"*: `1`
>    - *"Base name for new files in collection"*: `split_file`
>    - *"Method to allocate records to new files"*: `Maintain record order`
>
{: .hands_on}

> <hands-on-title>  </hands-on-title>
>
> 1. {% tool [Parse parameter value](param_value_from_file) %} with the following parameters:
>    - {% icon param-file %} *"Input file containing parameter to parse out of"*: 
>    - *"Select type of parameter to parse"*: `Text`
>    - *"Remove newlines ?"*: {% icon toggle %}  `Yes`
>
{: .hands_on}




> <hands-on-title> Cutting out SMILES column </hands-on-title>
>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/9.3+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: the tabular file in your history with the compound name and SMILES
>    - *"Operation"*: `Keep`
>    - *"Cut by"*: `fields`
>    - *"Delimited by"*: `Tab`
>    - *"Is there a header for the data’s columns"*: `No`
>    - *"List of Fields"*: `2`
>
{: .hands_on}


> <hands-on-title>  </hands-on-title>
>
> 1. {% tool [Split file to dataset collection](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_to_collection/split_file_to_collection/0.5.2) %} with the following parameters:
>    - {% icon param-file %} *"Select the file type to split"*: `Tabular`
>    - *"Tabular file to split"*: ``
>    - *"Number of header lines to transfer to new files"*: `0`
>    - *"Split by row or by a column?"*: `By row`
>    - *"Specify number of output files or number of records per file?"*: `Number of records per file (‘chunk mode’)`
>    - *"Chunk size"*: `1`
>    - *"Base name for new files in collection"*: `split_file`
>    - *"Method to allocate records to new files"*: `Maintain record order`
>
{: .hands_on}


> <hands-on-title>  </hands-on-title>
>
> 1. {% tool [Compound conversion](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_compound_convert/openbabel_compound_convert/3.1.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: split from SMILES
>    - *"Append the specified text after each molecule title"*: parse parameter value from name column
>
{: .hands_on}


> <hands-on-title>  </hands-on-title>
>
> 1. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/9.3+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Datasets to concatenate"*: output from the previous step
>
{: .hands_on}

# 3D Conformer generation & optimization

## Generate conformers
This step involves generating three-dimensional (3D) conformers for each molecule from the input SDF (Structure Data File). Conformers are different spatial arrangements of a molecule that result from rotations around single bonds. The number of conformers to generate can be specified as an input parameter, with a default value of 1 if not provided. This process is crucial for exploring the possible shapes and energies that a molecule can adopt. The output of this step is a file containing the generated 3D conformers.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Generate conformers](toolshed.g2.bx.psu.edu/repos/bgruening/ctb_im_conformers/ctb_im_conformers/1.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output` (Input dataset)
>    - *"Number of conformers to generate"*: `{'id': 1, 'output_name': 'output'}`
>
{: .hands_on}

## Molecular Format Conversion
Converts the generated conformers from the SDF format to Cartesian coordinate (XYZ) format. The XYZ format lists the atoms in a molecule and their respective 3D coordinates, which is a common format used in computational chemistry for further processing and analysis.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Compound conversion](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_compound_convert/openbabel_compound_convert/3.1.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: `outfile` (output of **Generate conformers** {% icon tool %})
>    - *"Output format"*: `XYZ cartesian coordinates format`
>    - *"Split multi-molecule files into a collection"*: {% icon toggle %} `Yes`
>    - *"Add hydrogens appropriate for pH"*: `7.0`
>
{: .hands_on}


## Molecular optimization

Performs semi-empirical optimization on the molecules using the Extended Tight-Binding (xTB) method. This step optimizes the geometry of the molecules to find the lowest energy conformation. The level of optimization accuracy to be used can be specified as an input parameter, *"Optimization Levels"*. The default quantum chemical method is GFN2-xTB.


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [xtb molecular optimization](toolshed.g2.bx.psu.edu/repos/recetox/xtb_molecular_optimization/xtb_molecular_optimization/6.6.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Atomic coordinates file"*: `file_outputs` (output of **Compound conversion** {% icon tool %})
>    - *"Optimization Levels"*: ``
>
{: .hands_on}

# QCxMS Spectra Prediction 

## Sub-step with **QCxMS neutral run**

Prepares the necessary input files for the QCxMS production runs. These files are required for running the QCxMS simulations, which will predict the mass spectra of the molecules. This step typically formats the optimized molecular data into a format that can be used for the production simulations. QC Method: Specifies the quantum chemical method to use (string, options: GFN1-xTB or GFN2-xTB).
Outputs:
.in output: Input file for the QCxMS production run (File).
.start output: Start file for the QCxMS production run (File).
.xyz output: Cartesian coordinate file for the QCxMS production run (File).


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [QCxMS neutral run](toolshed.g2.bx.psu.edu/repos/recetox/qcxms_neutral_run/qcxms_neutral_run/5.2.1+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Molecule 3D structure [.xyz]"*: `output` (output of **xtb molecular optimization** {% icon tool %})
>    - *"QC Method"*: ``
>
{: .hands_on}


## QCxMS production run

Calculates the mass spectra for each molecule using QCxMS (Quantum Chemistry and Mass Spectrometry). This simulation generates .res files, which contain the raw results of the mass spectra calculations. These results are essential for predicting how the molecules will appear in mass spectrometry experiments.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [QCxMS production run](toolshed.g2.bx.psu.edu/repos/recetox/qcxms_production_run/qcxms_production_run/5.2.1+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"in files [.in]"*: `coords1` (output of **QCxMS neutral run** {% icon tool %})
>    - {% icon param-file %} *"start files [.start]"*: `coords2` (output of **QCxMS neutral run** {% icon tool %})
>    - {% icon param-file %} *"xyz files [.xyz]"*: `coords3` (output of **QCxMS neutral run** {% icon tool %})
>
{: .hands_on}


## Filter failed datasets

Filters out any failed runs from the dataset to ensure only successful results are processed further. This step is important to maintain the integrity and quality of the data being analyzed in subsequent steps. The output is a file containing only the successful mass spectra results (File).


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter failed datasets](__FILTER_FAILED_DATASETS__) %} with the following parameters:
>    - {% icon param-file %} *"Input Collection"*: `res_files` (output of **QCxMS production run** {% icon tool %})
>
{: .hands_on}


## QCxMS get results

Converts the filtered .res files from the QCxMS production run into simulated mass spectra in MSP (Mass Spectrum Peak) file format. The MSP format is widely used for storing and sharing mass spectrometry data, enabling easy comparison and analysis of the results.


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [QCxMS get results](toolshed.g2.bx.psu.edu/repos/recetox/qcxms_getres/qcxms_getres/5.2.1+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Molecule 3D structure [.xyz]"*: `file_outputs` (output of **Compound conversion** {% icon tool %})
>    - {% icon param-file %} *"res files [.res]"*: `output` (output of **Filter failed datasets** {% icon tool %})
>
{: .hands_on}


# Conclusion
[key history](https://usegalaxy.eu/u/hechth/h/end-to-end-ei-mass-spectra-prediction-workflow-using-qcxms-1)
