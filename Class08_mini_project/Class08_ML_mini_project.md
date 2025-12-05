# Class 8: ML Mini Project (Breast Cancer Analysis)
Wade Ingersoll (PID: A69038080)

- [Background](#background)
- [Data import](#data-import)
- [Exploratory data analysis](#exploratory-data-analysis)

## Background

The goal of today’s mini-project is to explore a complete analysis using
the unsupervised learning techniques covered in the last class. We will
extend what we learned by combining PCA as a preprocessing step to
clustering using data that consist of measurements of cell nuclei of
human breast massess.

The data itself comes from the Wisconsin Breast Cancer Diagnostic Data
Set first reported by K. P. Benne and O. L. Mangasarian: “Robust Linear
Programming Discrimination of Two Linearly Inseparable Sets”.

Values in this data set describe characteristics of the cell nuclei
present in digitized images of a fine needle aspiration (FNA) of a
breast mass.

## Data import

The data is available as a CSV from class website:

``` r
wisc.df <- read.csv("WisconsinCancer (2).csv", row.names=1)
```

Make sure we do not include sample id or diagnosis columns in the data
that we analyze below.

``` r
diagnosis <- as.factor(wisc.df$diagnosis)
wisc.data <- wisc.df[, -1]
dim(wisc.data)
```

    [1] 569  30

## Exploratory data analysis

> Q1. How many observations are in this dataset?

There are 569 observations/samples/patients in the data set. I used
nrow(wisc.data) to find that number.

> Q2. How many of the observations have a malignant diagnosis?

There are 212; see code below (two ways to determine)

``` r
sum(wisc.df$diagnosis == "M")
```

    [1] 212

``` r
table(wisc.df$diagnosis)
```


      B   M 
    357 212 

> Q3. How many variables/features in the data are suffixed with \_mean?

There are 10; see code below

``` r
#colnames(wisc.data)
length( grep("_mean", colnames(wisc.data)) )
```

    [1] 10

Easier to read:

``` r
n <- colnames(wisc.data)
inds <- grep("_mean", n)
length(inds)
```

    [1] 10
