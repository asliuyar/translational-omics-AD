---
title: 'Synapse and AD Knowledge Portal'
teaching: 40
exercises: 10
---

:::::::::::::::::::::::::::::::::::::: questions 

- How to work with Synapse R client?
- How to work with data in AD Knowledge Portal?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Explain how to use Synapser Package.
- Demonstrate how to locate data and metadata in the Portal.
- Demonstrate how to download data from the Portal programmatically.

::::::::::::::::::::::::::::::::::::::::::::::::



## Working with AD Portal metadata 

We have now downloaded several metadata files and an RNAseq counts file from the
portal. For our next exercises, we want to read those files in as R data so we 
can work with them.

We can see from the `download_table` we got during the bulk download step
that we have five metadata files. Two of these should be the individual
and biospecimen files, and three of them are assay metadata files.


``` r
download_table %>%
  dplyr::select(name, metadataType, assay)
```

We are only interested in RNAseq data, so we will only read in the
individual, biospecimen, and RNAseq assay metadata files.


``` r
# counts matrix
counts <- read_tsv("data/htseqcounts_5XFAD.txt", 
                   show_col_types = FALSE)

# individual metadata
ind_meta <- read_csv("data/Jax.IU.Pitt_5XFAD_individual_metadata.csv", 
                     show_col_types = FALSE)

# biospecimen metadata
bio_meta <- read_csv("data/Jax.IU.Pitt_5XFAD_biospecimen_metadata.csv", 
                     show_col_types = FALSE)

# assay metadata
rna_meta <- read_csv("data/Jax.IU.Pitt_5XFAD_assay_RNAseq_metadata.csv", 
                     show_col_types = FALSE)
```

Let’s examine the data and metadata files a bit before we begin our
analyses. We start by exploring the `counts` data that we read in using the
tidyverse `read_tsv()` (*read*_*t*ab-*s*eparated *v*alues) function. This 
function reads data in as a *tibble*, a kind of data table with some nice 
features that avoid some bad habits of the base R `read.csv()` function. Calling 
a `tibble` object will print the first ten rows in a nice tidy output. Doing the 
same for a base R dataframe read in with `read.csv()` will print the whole thing 
until it runs out of memory. If you want to inspect a large dataframe, use 
`head(df)` to view the first several rows only.


``` r
counts
```

``` output
# A tibble: 55,489 × 73
   gene_id `32043rh` `32044rh` `32046rh` `32047rh` `32048rh` `32049rh` `32050rh`
   <chr>       <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
 1 ENSG00…     22554         0         0         0     16700         0         0
 2 ENSG00…    344489         4         0         1    260935         6         8
 3 ENSMUS…      5061      3483      3941      3088      2756      3067      2711
 4 ENSMUS…         0         0         0         0         0         0         0
 5 ENSMUS…       208       162       138       127        95       154       165
 6 ENSMUS…        44        17        14        28        23        24        14
 7 ENSMUS…       143        88       121       117       115       109        75
 8 ENSMUS…        22         6        10        11        11        19        24
 9 ENSMUS…      7165      5013      5581      4011      4104      5254      4345
10 ENSMUS…      3728      2316      2238      1965      1822      1999      1809
# ℹ 55,479 more rows
# ℹ 65 more variables: `32052rh` <dbl>, `32053rh` <dbl>, `32057rh` <dbl>,
#   `32059rh` <dbl>, `32061rh` <dbl>, `32062rh` <dbl>, `32065rh` <dbl>,
#   `32067rh` <dbl>, `32068rh` <dbl>, `32070rh` <dbl>, `32073rh` <dbl>,
#   `32074rh` <dbl>, `32075rh` <dbl>, `32078rh` <dbl>, `32081rh` <dbl>,
#   `32088rh` <dbl>, `32640rh` <dbl>, `46105rh` <dbl>, `46106rh` <dbl>,
#   `46107rh` <dbl>, `46108rh` <dbl>, `46109rh` <dbl>, `46110rh` <dbl>, …
```

The data file has a column of ENSEMBL `gene_id`s and then a bunch of columns 
with count data, where the column headers correspond to the `specimenID`s. These 
`specimenID`s should all be in the RNAseq assay metadata file, so let’s check.


``` r
# what does the RNAseq assay metadata look like?
rna_meta
```

``` output
# A tibble: 72 × 12
   specimenID platform   RIN   rnaBatch libraryBatch sequencingBatch libraryPrep
   <chr>      <chr>      <lgl>    <dbl>        <dbl>           <dbl> <chr>      
 1 32043rh    IlluminaN… NA           1            1               1 polyAselec…
 2 32044rh    IlluminaN… NA           1            1               1 polyAselec…
 3 32046rh    IlluminaN… NA           1            1               1 polyAselec…
 4 32047rh    IlluminaN… NA           1            1               1 polyAselec…
 5 32049rh    IlluminaN… NA           1            1               1 polyAselec…
 6 32057rh    IlluminaN… NA           1            1               1 polyAselec…
 7 32061rh    IlluminaN… NA           1            1               1 polyAselec…
 8 32065rh    IlluminaN… NA           1            1               1 polyAselec…
 9 32067rh    IlluminaN… NA           1            1               1 polyAselec…
10 32070rh    IlluminaN… NA           1            1               1 polyAselec…
# ℹ 62 more rows
# ℹ 5 more variables: libraryPreparationMethod <lgl>, isStranded <lgl>,
#   readStrandOrigin <lgl>, runType <chr>, readLength <dbl>
```


``` r
# are all the column headers from the counts matrix 
# (except the first "gene_id" column) in the assay metadata?
all(colnames(counts[-1]) %in% rna_meta$`specimenID`)
```

``` output
[1] TRUE
```

**Assay metadata** 

The assay metadata contains information about how data
was generated on each sample in the assay. Each `specimenID` represents a
unique sample. We can use some tools from `dplyr` to explore the metadata.


``` r
# how many unique specimens were sequenced?
n_distinct(rna_meta$`specimenID`)
```

``` output
[1] 72
```


``` r
# were the samples all sequenced on the same platform?
distinct(rna_meta, platform)
```

``` output
# A tibble: 1 × 1
  platform           
  <chr>              
1 IlluminaNovaseq6000
```


``` r
# were there multiple sequencing batches reported?
distinct(rna_meta, sequencingBatch) 
```

``` output
# A tibble: 1 × 1
  sequencingBatch
            <dbl>
1               1
```

**Biospecimen metadata** 

The biospecimen metadata contains specimen-level information, including organ 
and tissue the specimen was taken from, how it was prepared, etc. Each 
`specimenID` is mapped to an `individualID`.


``` r
# all specimens from the RNAseq assay metadata file should be in the biospecimen file
all(rna_meta$`specimenID` %in% bio_meta$`specimenID`)
```

``` output
[1] TRUE
```


``` r
# but the biospecimen file also contains specimens from different assays
all(bio_meta$`specimenID` %in% rna_meta$`specimenID`)
```

``` output
[1] FALSE
```

**Individual metadata** 

The individual metadata contains information about all the individuals in the 
study, represented by unique `individualID`s. For humans, this includes 
information on age, sex, race, diagnosis, etc. For MODEL-AD mouse models, the
individual metadata has information on model genotypes, stock numbers, diet, and 
more.


``` r
# all `individualID`s in the biospecimen file should be in the individual file
all(bio_meta$`individualID` %in% ind_meta$`individualID`)
```

``` output
[1] TRUE
```


``` r
# which model genotypes are in this study?
distinct(ind_meta, genotype)
```

``` output
# A tibble: 2 × 1
  genotype        
  <chr>           
1 5XFAD_carrier   
2 5XFAD_noncarrier
```

**Joining metadata** 

We use the three-file structure for our metadata because it allows us to store 
metadata for each study in a tidy format. Every line in the assay and 
biospecimen files represents a unique specimen, and every line in the individual 
file represents a unique individual. This means the files can be easily joined 
by `specimenID` and `individualID` to get all levels of metadata that apply to a
particular data file. We will use the `left_join()` function from the `dplyr` 
package, and the `%\>%` operator from the `magrittr` package. If you are 
unfamiliar with the pipe, think of it as a shorthand for “take this (the 
preceding object) and do that (the subsequent command)”. See [magrittr](https://magrittr.tidyverse.org/) for more info on piping in R.


``` r
# join all the rows in the assay metadata
# that have a match in the biospecimen metadata
joined_meta <- rna_meta %>%               # start with the rnaseq assay metadata
  left_join(bio_meta, by = "specimenID") %>%  # join rows from biospecimen 
                                                # that match specimenID 
  left_join(ind_meta, by = "individualID") # join rows from individual 
                                              # that match individualID

joined_meta
```

``` output
# A tibble: 72 × 53
   specimenID platform   RIN   rnaBatch libraryBatch sequencingBatch libraryPrep
   <chr>      <chr>      <lgl>    <dbl>        <dbl>           <dbl> <chr>      
 1 32043rh    IlluminaN… NA           1            1               1 polyAselec…
 2 32044rh    IlluminaN… NA           1            1               1 polyAselec…
 3 32046rh    IlluminaN… NA           1            1               1 polyAselec…
 4 32047rh    IlluminaN… NA           1            1               1 polyAselec…
 5 32049rh    IlluminaN… NA           1            1               1 polyAselec…
 6 32057rh    IlluminaN… NA           1            1               1 polyAselec…
 7 32061rh    IlluminaN… NA           1            1               1 polyAselec…
 8 32065rh    IlluminaN… NA           1            1               1 polyAselec…
 9 32067rh    IlluminaN… NA           1            1               1 polyAselec…
10 32070rh    IlluminaN… NA           1            1               1 polyAselec…
# ℹ 62 more rows
# ℹ 46 more variables: libraryPreparationMethod <lgl>, isStranded <lgl>,
#   readStrandOrigin <lgl>, runType <chr>, readLength <dbl>,
#   individualID <dbl>, specimenIdSource <chr>, organ <chr>, tissue <chr>,
#   BrodmannArea <lgl>, sampleStatus <chr>, tissueWeight <lgl>,
#   tissueVolume <lgl>, nucleicAcidSource <lgl>, cellType <lgl>,
#   fastingState <lgl>, isPostMortem <lgl>, samplingAge <lgl>, …
```

We now have a very wide dataframe that contains all the available
metadata on each specimen in the RNAseq data from this study. This
procedure can be used to join the three types of metadata files for
every study in the AD Knowledge Portal, allowing you to filter
individuals and specimens as needed based on your analysis criteria!


``` r
library(lubridate)

# convert columns of strings to month-date-year format
joined_meta_time <- joined_meta %>% 
  mutate(dateBirth = mdy(dateBirth), dateDeath = mdy(dateDeath)) %>% 
  
  # create a new column that subtracts dateBirth from dateDeath in days, 
  # then divide by 30 to get months
  mutate(timepoint = as.numeric(difftime(dateDeath, dateBirth, units ="days"))/30) %>% 
  
  # convert numeric ages to timepoint categories
  mutate(timepoint = case_when(timepoint > 10 ~ "12 mo",
                               timepoint < 10 & timepoint > 5 ~ "6 mo",
                               timepoint < 5 ~ "4 mo"))

covars_5XFAD <- joined_meta_time %>%
  dplyr::select(individualID, specimenID, sex, genotype, timepoint) %>% 
  distinct() %>% 
  as.data.frame()
rownames(covars_5XFAD) <- covars_5XFAD$specimenID

head(covars_5XFAD)
```

``` output
        individualID specimenID    sex         genotype timepoint
32043rh        32043    32043rh female    5XFAD_carrier     12 mo
32044rh        32044    32044rh   male 5XFAD_noncarrier     12 mo
32046rh        32046    32046rh   male 5XFAD_noncarrier     12 mo
32047rh        32047    32047rh   male 5XFAD_noncarrier     12 mo
32049rh        32049    32049rh female 5XFAD_noncarrier     12 mo
32057rh        32057    32057rh female 5XFAD_noncarrier     12 mo
```

We will save `joined_meta` for the next lesson.


``` r
saveRDS(covars_5XFAD, file = "data/covars_5XFAD.rds")
```

## Single Specimen files

For files that contain data from a single specimen (e.g. raw sequencing files, 
raw mass spectra, etc.), we can use the Synapse annotations to associate these 
files with the appropriate metadata.

::::::::::::::::::::::::::::::::::::: challenge

## Exercise

Use *Explore Data* to find all RNAseq files from the `Jax.IU.Pitt_5XFAD` study. 
If we filter for data where `Study = "Jax.IU.Pitt_5XFAD"` and `Assay = "rnaSeq"`
we will get a list of 148 files, including raw fastqs and processed counts data.

:::::::::::::::: solution

:::::::::::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::::::::


We can use the function `synGetAnnotations` to view the annotations associated 
with any file without downloading the file.


``` r
# the synID of a random fastq file from this list
random_fastq <- "syn22108503"

# extract the annotations as a nested list
fastq_annotations <- synGetAnnotations(random_fastq)

fastq_annotations
```

The file annotations let us see which study the file is associated with (`Jax.IU.Pitt.5XFAD`), which species it is from (mouse), which assay generated 
the file (rnaSeq), and a whole bunch of other properties. Most importantly, single-specimen files are annotated with the `specimenID` of the specimen in the 
file, and the `individualID` of the individual that specimen was taken from. We 
can use these annotations to link files to the rest of the metadata, including 
metadata that is not in annotations. This is especially helpful for human 
studies, as potentially identifying information like age, race, and diagnosis is 
not included in file annotations.


``` r
# find records belonging to the individual this file maps
# to in our joined metadata
joined_meta %>% 
  filter(`individualID` == fastq_annotations$individualID[[1]])
```

## Annotations during bulk download

When bulk downloading many files, the best practice is to preserve the download manifest that is generated which lists all the files, their `synIDs`, and all 
their annotations. If using the Synapse R client, follow the instructions in the 
Bulk download files section above.

If we use the “Programmatic Options” tab in the AD Portal download menu to 
download all 148 rnaSeq files from the 5XFAD study, we would get a table query 
that looks like this:


``` r
query <- synTableQuery("SELECT * FROM syn11346063.37 
                       WHERE ( ( \"study\" HAS ( 'Jax.IU.Pitt_5XFAD' ) ) 
                       AND ( \"assay\" HAS ( 'rnaSeq' ) ) )")
```

As we saw previously, this downloads a csv file with the results of our AD 
Portal query. Opening that file lets us see which specimens are associated with
which files:


``` r
annotations_table <- read_csv(query$filepath, show_col_types = FALSE)

annotations_table
```

You could then use 
`purrr::walk(download_table$id, ~synGet(.x, downloadLocation = ))` to walk 
through the column of `synIDs` and download all 148 files. However, because 
these are large files, it might be preferable to use the Python client or 
command line client for increased speed.

Once you’ve downloaded all the files in the `id` column, you can link those 
files to their annotations by the `name` column.

We'll use the "random fastq" that we got annotations for earlier to avoid 
downloading the whole 3GB file, we'll use synGet with `downloadFile = FALSE` to 
get only the Synapse entity information, rather than the file. If we downloaded 
the actual file, we could find it in the directory and search using the 
filename. Since we're just downloading the Synapse entity wrapper object, we'll 
use the file name listed in the object properties.


``` r
fastq <- synGet(random_fastq, downloadFile = FALSE)

# filter the annotations table to rows that match the fastq filename
annotations_table %>% 
  filter(name == fastq$properties$name)
```

## Multispecimen files

Multispecimen files in the AD Knowledge Portal are files that contain data or information from more than one specimen. They are not annotated with 
`individualID`s or `specimenID`s, since these files may contain numbers of 
specimens that exceed the annotation limits. These files are usually processed 
or summary data (gene counts, peptide quantifications, etc), and are always 
annotated with `isMultiSpecimen = TRUE`.

If we look at the processed data files in the table of 5XFAD RNAseq file 
annotations we just downloaded, we will see that it `isMultiSpecimen = TRUE`, 
but `individualID` and `specimenID` are blank:


``` r
annotations_table %>% 
  filter(fileFormat == "txt") %>% 
  dplyr::select(name, individualID, specimenID, isMultiSpecimen)
```

The multispecimen file should contain a row or column of `specimenID`s that 
correspond to the `specimenID`s used in a study’s metadata, as we have seen with 
the 5XFAD counts file.

In this example, we take a slice of the counts data to reduce computation, 
transpose it so that each row represents a single specimen, and then join it 
to the joined metadata by the `specimenID`.


``` r
counts %>% 
  slice_head(n = 5) %>% 
  t() %>% 
  as_tibble(rownames = "specimenID") %>% 
  left_join(joined_meta, by = "specimenID")
```

``` output
# A tibble: 73 × 58
   specimenID V1    V2    V3    V4    V5    platform RIN   rnaBatch libraryBatch
   <chr>      <chr> <chr> <chr> <chr> <chr> <chr>    <lgl>    <dbl>        <dbl>
 1 gene_id    "ENS… "ENS… "ENS… "ENS… "ENS… <NA>     NA          NA           NA
 2 32043rh    " 22… "344… "  5… "   … "   … Illumin… NA           1            1
 3 32044rh    "   … "   … "348… "   … " 16… Illumin… NA           1            1
 4 32046rh    "   … "   … "394… "   … " 13… Illumin… NA           1            1
 5 32047rh    "   … "   … "308… "   … " 12… Illumin… NA           1            1
 6 32048rh    " 16… "260… "  2… "   … "   … Illumin… NA           1            1
 7 32049rh    "   … "   … "306… "   … " 15… Illumin… NA           1            1
 8 32050rh    "   … "   … "271… "   … " 16… Illumin… NA           1            1
 9 32052rh    " 19… "337… "  3… "   … "   … Illumin… NA           1            1
10 32053rh    " 14… "206… "  3… "   … "   … Illumin… NA           1            1
# ℹ 63 more rows
# ℹ 48 more variables: sequencingBatch <dbl>, libraryPrep <chr>,
#   libraryPreparationMethod <lgl>, isStranded <lgl>, readStrandOrigin <lgl>,
#   runType <chr>, readLength <dbl>, individualID <dbl>,
#   specimenIdSource <chr>, organ <chr>, tissue <chr>, BrodmannArea <lgl>,
#   sampleStatus <chr>, tissueWeight <lgl>, tissueVolume <lgl>,
#   nucleicAcidSource <lgl>, cellType <lgl>, fastingState <lgl>, …
```


::::::::::::::::::::::::::::::::::::: keypoints 

- 

::::::::::::::::::::::::::::::::::::::::::::::::

