# HLA-ClassI-Tumorigenesis
![](https://i.imgur.com/552qtze.png)

 HLA-I is an important molecule in killing tumor cells along with CD8+ T cells. CD8+ T cells recognize and kill tumor cells with HLA-I tumor antigens in early tumorigenesis, the efficiency of which differs according to antigen recognition coverage. We assessed the tumor type- and driver mutation-specificity in the association between tumor onset age and HLA-I zygosity. 

<br/>

## Dataset
TCGA and ICGC data were used for pan-oncological analysis.

<br/>

### TCGA
The Cancer Genome Atlas (TCGA) is a dataset that includes large-scale genome sequencing. Some data in TCGA is controlled data and requires dbGaP access `(dbGaP accession number: phs000178)`

<br/>

| Contents               | Source                                                                                                                                                       | File                                                                                       |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------ |
| Clinical data          | [pancanatlas](https://gdc.cancer.gov/about-data/publications/pancanatlas)                                                                                    | TCGA-Clinical Data Resource (CDR) Outcome                                                  |
| HLA                    | [optitype](https://gdc.cancer.gov/about-data/publications/panimmune), [Polysolver](https://gdc.cancer.gov/about-data/publications/PanCanAtlas-Splicing-2018) | HLA class I type calls from tumor RNA-seq using fastq as input and HLA types, respectively |
| Mutation               | [Variant Calls](https://gdc.cancer.gov/about-data/publications/mc3-2017)                                                                                     | MC3 Controlled MAF                                                                         |
| Copy number alteration | [Gistic2](https://gdac.broadinstitute.org/runs/analyses__2016_01_28/reports/cancer/KIRC-TP/CopyNumber_Gistic2/nozzle.html)                                   | All Data By Genes                                                                          |
| Viral status           | [panimmune](https://gdc.cancer.gov/about-data/publications/panimmune)                                                                                        | Viral Read Counts                                                                          |
| Pathogen               | [CharGer](https://gdc.cancer.gov/about-data/publications/PanCanAtlas-Germline-AWG)                                                                           | Prioritized, cancer related variants discovered in 10,389 cases                            |
| RNA | [Broad GDAC Firehose](http://firebrowse.org/?cohort=KIRC&download_dialog=true) | mRNAseq_Preprocess |

---
<br/>

###  ICGC dataset
The International Cancer Genome Consortium (ICGC) provides genomic raw and analyze data from 50 different cancer types. Controlled data can be requested through The Data Access Compliance Office (DACO).

To Validate the effect of VHL loss, a RECA-EU cohort with histology type(`8310/3`) was used.

<br/>

| Contents     | Source                                                      | File                      |
| ------------ | ----------------------------------------------------------- | ------------------------- |
| Raw data     | [ICGC](https://dcc.icgc.org/)                               |                           |
| Allelic loss | [DCC](https://dcc.icgc.org/releases/PCAWG/driver_mutations) | panorama driver mutations |
| HLA  | [DCC](https://dcc.icgc.org/releases/PCAWG/hla_and_neoantigen) | 4 digit normal v2 |
---
