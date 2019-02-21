class14
================

``` r
read.sample <- read.csv("373510-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head(read.sample)
```

    ##   Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s.
    ## 1                  HG00551 (F)                       G|A ALL, AMR, PUR
    ## 2                  HG00553 (M)                       A|A ALL, AMR, PUR
    ## 3                  HG00554 (F)                       A|A ALL, AMR, PUR
    ## 4                  HG00637 (M)                       G|G ALL, AMR, PUR
    ## 5                  HG00638 (F)                       G|A ALL, AMR, PUR
    ## 6                  HG00640 (M)                       A|G ALL, AMR, PUR
    ##   Father Mother
    ## 1      -      -
    ## 2      -      -
    ## 3      -      -
    ## 4      -      -
    ## 5      -      -
    ## 6      -      -

``` r
table(read.sample$Genotype..forward.strand.)
```

    ## 
    ## A|A A|G G|A G|G 
    ## 122  85  80  60

``` r
sum(read.sample$Genotype..forward.strand.== "G|G")
```

    ## [1] 60

``` r
(table(read.sample$Genotype..forward.strand.)/nrow(read.sample))*100
```

    ## 
    ##      A|A      A|G      G|A      G|G 
    ## 35.15850 24.49568 23.05476 17.29107

``` r
q <- "DDDDCDEDCDDDDBBDDDCC@"
```

``` r
q <- "DDDDCDEDCDDDDBBDDDCC@"
library(seqinr)
library(gtools)
asc(s2c(q))-33
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

``` r
geno <- read.table("rs8067378_ENSG00000172057.6.txt")
head(geno)
```

    ##    sample geno      exp
    ## 1 HG00367  A/G 28.96038
    ## 2 NA20768  A/G 20.24449
    ## 3 HG00361  A/A 31.32628
    ## 4 HG00135  A/A 34.11169
    ## 5 NA18870  G/G 18.25141
    ## 6 NA11993  A/A 32.89721

``` r
table(geno$geno)
```

    ## 
    ## A/A A/G G/G 
    ## 108 233 121

``` r
inds <- geno$geno == "G/G"
summary(geno$exp[inds])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   6.675  16.903  20.074  20.594  24.457  33.956

``` r
inds <- geno$geno == "A/G"
summary(geno$exp[inds])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   7.075  20.626  25.065  25.397  30.552  48.034

``` r
inds <- geno$geno == "G/A"
summary(geno$exp[inds])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 

``` r
inds <- geno$geno == "A/A"
summary(geno$exp[inds])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   11.40   27.02   31.25   31.82   35.92   51.52

``` r
boxplot(exp~geno, geno)
```

![](class14_files/figure-markdown_github/unnamed-chunk-10-1.png)
