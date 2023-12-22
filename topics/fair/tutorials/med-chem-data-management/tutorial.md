---
layout: tutorial_hands_on
title: Data management in Medicinal Chemistry

zenodo_link: ''

questions:
- Why does the medicinal chemistry research produce so much data?
- How can big data be useful for medicinal chemists?
- What are the current publicly available databases that can be used for drug discovery and development?

objectives:
- Learn some terminology from the field of medicinal chemistry
- Understand the idea of data-driven medicinal chemistry
- Explore the databases used for drug discovery and development

time_estimation: 1H

key_points:
- The power of big data might be useful to shape innovations of the future in the pharmaceutical field.
- There are many publicly available databases used for drug discovery and development and they might look at the same medicines from different angles and classify them based on various factors.
  
tags:
- fair
- data management
- medicinal-chemistry
- computational-chemistry
  
priority: 6

contributions:
  authorship:
    - wee-snufkin

subtopic: fair-data

requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro
      - data-management
   
funding: 
  - elixir-fair-data

follow_up_training:
  - type: "internal"
    topic_name: computational-chemistry
---

# Introduction 
The development of medicinal chemistry is more and more advancing. Big pharmaceutical companies, research institutes and universities are working on ground-breaking solutions to help patients combat all kinds of diseases. During that development process, tons of data are generated – not only from the lab environment but also from clinical trials. Given that the discovery of more potent, safer and cheaper drugs is the ultimate goal of all research bodies, we should all focus on making the data we gather FAIR: Findable, Accessible, Interoperable, and Reusable to push the boundaries of drug development even further. 

With the currently available methods such as artificial intelligence, machine learning, many toolkits, software and access to various databases, it seems that managing big data is now inherently linked to medicinal chemistry and hence ensures that this area is as efficient as it can be. 


# Principles of medicinal chemistry research

## Drug likeliness

## Ten Vs of big data for drug discovery
Volume: size of data. 
Velocity: speed of new data generation. 
Variety: various formats of data. 
Veracity: quality of data. 
Validity: authenticity of data. 
Vocabulary: terminology of data. 
Venue: platform of data generation. 
Visualization: view of data. 
Volatility: duration of data usefulness. 
Value: potential of data usefulness to reduce the cost of drug discovery and development.


## The Lipinski rule of 5
There are many factors that medicinal chemists take into account while designing new drugs. 
The rule of thumb developed back in 1997 by Christopher Lipinski and colleagues tries to predict the likelihood that a given small molecule can be orally active. This Lipinski rule of 5 favours molecules as potential oral drug candidates if:
-	The molecular mass is less than 500
-	calculated logarithm of the octanol−water partition coefficient (clogP) is less than 5
-	they have up to 5 H-bond donors
-	they have up to 10 H-bond acceptors

Read more about re-assessing the rule of 5: {% cite Mullard2018 %}.

However, nowadays there are more and more drugs being developed which don’t comply with those rules and regardless are still effective. There are voices from the scientific community, pointing out that “We are in danger of repeating our past mistakes if we assume these new modalities are not 'drug-like' and cannot be oral drugs because they are not [rule of 5] compliant” (Michael Shultz from Novartis {% cite Mullard2018 %}). Then, in {% cite OHagan2014 %} we read "This famous "rule of 5" has been highly influential in this regard, but only about 50 % of orally administered new chemical entities actually obey it." 


## Veber's rule

In addition to the rule of 5, there are also other criteria, one of which is the Veber's rule, which helps to predict good oral bioavailability based on:
- 10 or fewer rotatable bonds
- Polar surface area no greater than 140 Å <sup>2</sup>


# FAIR MedChem
Those are the properties that the compounds can be easily searched by if all the information is included in the database. There is a need of assessing other molecular properties. This is one of the reasons why new repositories are being developed – they are more specific and gather particular properties of interest [EXAMPLE?]. 
How much easier the life of scientists could be if the relevant data was publicly available, well-ordered and contained the needed metadata? By submitting the data with necessary information to the repository is a good way to make the data FAIR. This will make it:
- **F**indable, as the data will be given specific identifiers
- **A**ccessible, as the data will be available online, open and free where possible
-	**I**nteroperable, as the repository will often enforce the use of formalised, consistent language
-	**R**eusable, as the data will be released under a license with detailed provenance
In the repositories 



# Chemical and pharmacological databases

Currently there are lots of publicly available databases storing information about hit molecules, chemical structures and properties, drug targets, pharmaceuticals… In the recent paper {% cite Zhao2020 %} the authors collated relevant databases to advance computer-aided drug discovery and divided them into several main groups:
- Chemical collections
- Drug / drug-like compounds
- Drug targets, including genomics and proteomics data
- Biological data from assay screening, metabolism, efficacy studies
- Drug liabilities and toxicities
- Clinical databases

Below you will find the databases listed under corresponding categoories, all taken from {% cite Zhao2020 %} paper. 

There are many more databases available, and many are still being developed. They are often quite specific and contain certain types of compounds (eg. PROTACS) or are aimed at particular disease. If the listed databases are not specific enough for your research, you can try luck by searching smaller, more specific databases. 

## Chemical collections

| Database | Description | Size (as of 29 October 2019) |
|:---:|:---:|:---:|
| [Enamine REAL Database](https://enamine.net/hit-finding/compound-collections/real-database) | Tool used to find new hit molecules using largescale virtual screening   and for searching analogs to hit molecules | >700 million compounds that comply with ‘rule of 5’ and Veber   criteria |
| [ZINC](http://zinc.docking.org/) | Contains compound information including 2D/3D structure, purchasability,   target, and biologyrelated information | >230 million compounds in 3D formats and >750 million compounds for   analogsearching |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov) | Contains chemical molecule (mostly small molecule) information, including   chemical structures, identifiers, chemical and physical properties,   biological activities, safety and toxicity data | 97 million compounds, 236 million substances, 268 million bioactivities |
| [ChemSpider](www.chemspider.com) | Free chemical structure database providing fast access to >78 million   structures, properties, and associated information | >78 million compound structures |
| [SCUBIDOO](http://kolblab.org/scubidoo/index.php) | Freely accessible database that currently holds 21 million virtual   products originating from small library of building blocks and collection of   robust organic reactions | 21 million virtual products |
| [ChEMBL](www.ebi.ac.uk/chembl) | Manually curated database of bioactive molecules with drug-like   properties; brings together chemical, bioactivity, and genomic data to aid   translation of genomic information into effective new drugs. | >1.9 million compounds, 1.1 million pieces of assay information |
| [TCM-Mesh](http://mesh.tcm.microbioinformatics.org/) | Integration of a database and a data-mining system for network   pharmacology analysis of all respects of traditional Chinese medicine,   including herbs, herbal ingredients, targets, related diseases, adverse   effect, and toxicity | 383 840 compounds, 6235 herbs |
| [Super Natural II](http://bioinf-applied.charite.de/supernatural_new/index.php) | Contains natural compounds, including information about corresponding 2D   structures, physicochemical properties, predicted toxicity class and   potential vendors | 325 508 natural compounds |
| [BIAdb](https://webs.iiitd.edu.in/raghava/biadb/index.html) | Comprehensive database of benzylisoquinoline alkaloids, containing   information about 846 unique benzylisoquinoline alkaloids | ~846 unique benzylisoquinoline alkaloids |

## Drug / drug-like compounds

| Database | Description | Size (as of 29 October 2019) |
|:---:|:---:|:---:|
| [AICD](http://956023.ichengyun.net/AICD/index.php) | Anti-Inflammatory Compounds Database (AICD) deposits compounds with   potential antiinflammation activities | 79 781 small molecules |
| [DrugBank](www.drugbank.ca) | Unique bioinformatics and cheminformatics resource that combines detailed   drug data with comprehensive drug target information | 13 441 drug entries |
| [ReFRAME](https://reframedb.org/) | Screening library of 12 000 molecules assembled by combining three   databases (Clarivate Integrity, GVK Excelra GoStar, and Citeline   Pharmaprojects) to facilitate drug repurposing | 12 000 molecules |
| [SuperDRUG2](http://cheminfo.charite.de/superdrug2/) | Contains approved/marketed drugs with regulatory details, chemical   structures (2D and 3D), dosage, biological targets, physicochemical   properties, external identifiers, adverse effects, and PK data | >4600 active pharmaceutical ingredients |
| [Drugs@FDA database](www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files) | Information about drugs from FDA | ~23 391 drug application records |
| [e-Drug3D](https://chemoinfo.ipmc.cnrs.fr/MOLDB/index.php) | Contains 1930 molecular structures approved by FDA between 1939 and 2019   with a molecular weight <2000 | 1930 drugs |

## Drug targets, including genomics and proteomics data

| Database | Description | Size (as of 29 October 2019) |
|:---:|:---:|:---:|
| [BindingDB](www.bindingdb.org/bind/index.jsp) | Public, web-accessible database of   measured binding affinities, focusing chiefly on interactions of proteins   considered to be candidate drug targets with ligands that are small,   drug-like molecules |  |
| [Supertarget](http://insilico.charite.de/supertarget/index.php?site=home) | An extensive web resource for   analyzing drugtarget interactions. |  |
| [Ligand Expo](http://ligand-expo.rutgers.edu/index.html) | Provides chemical and structural   information about small molecules within structure entries of Protein Data   Bank. |  |
| [PDBeChem](www.ebi.ac.uk/pdbe-srv/pdbechem/) | Consistent and enriched library of   ligands, small molecules, and monomers referenced as residues and hetgroups   in PDB entries |  |
| [PDBbind-CN](www.pdbbind-cn.org/) | Provides essential linkage between   energetic and structural information of biomolecular complexes, which is   helpful for various computational and statistical studies on molecular   recognition in biological systems |  |
| [STITCH](http://stitch.embl.de/) | Database integrating information   about interactions from metabolic pathways, crystal structures, binding   experiments, and drug–target relationships |  |
| [BioGRID](https://thebiogrid.org/) | The Biological General Repository   for Interaction Datasets is an open-access database on protein, genetic, and   chemical interactions for humans and all major model organisms |  |
| [Binding MOAD](http://bindingmoad.org/) | Created from a subset of Protein   Data Bank (PDB), containing every high-quality example of ligandprotein   binding. |  |
| [GPCRdb](www.gpcrdb.org) | Contains data from GPCRs, including   crystal structures, sequence alignments, and receptor mutations; can be   visualized in interactive diagrams; provides online analysis tools |  |
| [Guide to Pharmacology](www.guidetopharmacology.org/) | IUPHAR/BPS guide to pharmacology is   an openaccess, expert-curated database of molecular interactions between   ligands and their targets. |  |
| [GLASS](https://zhanglab.ccmb.med.umich.edu/GLASS/) | GPCR-Ligand Association (GLASS)   database is a manually curated repository for experimentally validated   GPCR–ligand interactions; along with relevant GPCR and chemical information,   GPCRligand association data are extracted and integrated into GLASS from literature   and public databases |  |

## Biological data from assay screening, metabolism, efficacy studies

| Database | Description | Size (as of 29 October 2019) |
|:---:|:---:|:---:|
| [HMDB](www.hmdb.ca/about) | Freely available electronic database containing detailed information   about small-molecule metabolites found in human body | 114 162 metabolite entries |
| [SMPDB](http://smpdb.ca/) | Small Molecule Pathway Database (SMPDB) is an interactive, visual   database containing >30 000 small-molecule pathways found in humans only | >30 000 small-molecule pathways |
| [TTD](http://db.idrblab.net/ttd/) | Therapeutic Target Database (TTD) is a database providing information   about known and explored therapeutic protein and nucleic acid targets,   targeted disease, pathway information and corresponding drugs directed at   each of these targets | 2589 targets, and 31 614 drugs |
| [BioCyc](https://biocyc.org/) | Collection of 7615 pathway/genome databases; each database in BioCyc   collection describes genome and metabolic pathways of a single organism | 7615 pathway/genome databases |
| [BiGG](http://bigg.ucsd.edu/) | Metabolic reconstruction of human metabolism designed for systems biology   simulation and metabolic flux balance modeling | 2004 proteins, 2766 metabolites, and 3311 metabolic and transport   reactions |
| [BRENDA](www.brenda-enzymes.org/) | Main collection of enzyme functional data available to scientific   community | At least 40 000 different enzymes from >6900 different organisms |
| [Reactome](https://reactome.org/) | Curated, peer-reviewed knowledgebase of biological pathways, including   metabolic pathways, and protein trafficking and signaling pathways | >9600 proteins, 9800 reactions, and 2000 pathways for humans |
| [BioModels Database](www.ebi.ac.uk/biomodels-main/) | Repository of computational models of biological processes; models   described from literature are manually curated and enriched with   crossreferences | 6753 patient-derived genomescale metabolic models, 112 898 metabolic   models etc. |
| [KEGG](https://www.genome.jp/kegg) | Database resource that integrates genomic, chemical, and systemic   functional information | 18 652 metabolites |
| [CARLSBAD](http://carlsbad.health.unm.edu/carlsbad/?mode=home) | Database and knowledge inference system that integrates multiple   bioactivity data sets to provide researchers with novel capabilities for   mining and exploration of available structure activity relationships (SAR)   throughout chemical biology space. | 932 852 CARLSBAD activities, 890 323 unique structure–target pairs, 3542   targets, 435 343 unique structures |
| [WOMBAT](http://dud.docking.org/wombat/) | Contains 331 872 entries, representing 1966 unique targets, with   bioactivity annotations | 268 246 unique structures |
| [Open NCI Database](https://cactus.nci.nih.gov/ncidb2.2/) | Maintained by the National Cancer Institute; contains small-molecule   information such as names, biological activities, structures; useful resource   for researchers working in cancer/AIDS fields | >250 000 compounds |
| [NPACT](http://crdd.osdd.net/raghava/npact/) | Provides information on plant-derived natural compound, including   structure, properties (physical, elemental, and topological), cancer type,   cell lines, inhibitory values (IC50, ED50, EC50, GI50), molecular targets,   commercial suppliers, and drug likeness of compounds | 1574 entries |
| [PKPB_DB](https://cfpub.epa.gov/ncea/risk/recordisplay.cfm?deid=204443) | Contains physiological parameter values for humans from early childhood   through senescence; intended to be used in physiologically based (PB) PK   modeling; also contains similar data for animals (primarily rodents) | N/A |

## Drug liabilities and toxicities

| Database | Description | Size (as of 29 October 2019) |
|:---:|:---:|:---:|
| [T3DB](www.t3db.ca/) | Unique bioinformatics resource that combines detailed toxin data with   comprehensive toxin target information | 3678 toxins |
| [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix/) | One of world’s largest toxicogenomic reference resources | ~ 600 drug molecules and 10 000 genes |
| [ACToR](https://actor.epa.gov/actor/home.xhtml) | Includes computational toxicology information about compounds, including   HTS, chemical exposure, sustainable chemistry (chemical structures and   physicochemical properties) and virtual tissue data | >500 000 chemicals |
| [SkinSensDB](https://cwtung.kmu.edu.tw/skinsensdb/) | Contains curated data from published AOP-related skin sensitization   assays | 710 unique chemicals |
| [SIDER](http://sideeffects.embl.de/download/) | Contains information on marketed medicines and their recorded adverse   drug reactions, including frequency, drug and adverse effect classifications | 1430 drugs with 5868 side effect information |
| [LTKB Benchmark Dataset](www.fda.gov/science-research/) | Contains drugs with potential to cause druginduced liver injury in   humans; established using FDA-approved prescription drug labels;   liver-toxicity-knowledge-base-ltkb/ ltkb-benchmark-dataset | 287 prescription drugs |
| [CTD](http://ctdbase.org/) | Comparative Toxicogenomics Database (CTD) is a premier public resource   for literature-based, manually curated associations between chemicals, gene   products, phenotypes, diseases, and environmental exposures | 13 378 unique chemicals and related information |

## Clinical databases

| Database | Description | Size (as of 29 October 2019) |
|:---:|:---:|:---:|
| [ClinicalTrials.gov](https://clinicaltrials.gov/) | Database of privately and publicly funded clinical studies conducted   around the world | ~ 324 429 research studies in all 50 US states and 209 countries |
| [AACT database](https://aact.ctti-clinicaltrials.org/) | Publicly available relational database that contains all information   (protocol and result data elements) about every study registered in   ClinicalTrials.gov. Content is downloaded from ClinicalTrials.gov daily and   loaded into AACT | ~ 324 429 research studies in all 50 US states and 209 countries |
| [EORTC Clinical Trials Database](www.eortc.org/clinical-trials/) | Contains information about EORTC clinical trials and clinical trials from   other organizations with EORTC participation | N/A |
| [Exposome-Explorer](http://exposome-explorer.iarc.fr/) | Contains detailed information on nature of biomarkers, populations and   subjects where measured, samples analyzed, methods used for biomarker   analyses, concentrations in biospecimens, correlations with external exposure   measurements, and biological reproducibility over time | 908 biomarkers |
| [PharmaGKB](www.pharmgkb.org/) | A pharmacogenomics knowledge resource that encompasses clinical   information about drug molecules | 733 drugs with their clinical information |

# Data management in chemistry using Galaxy

# Looking into the future: data-driven medicinal chemistry
