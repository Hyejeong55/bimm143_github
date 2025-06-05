# Class 7: Machine Learning 1
Hyejeong Choi (PID: A16837133)

- [Clustering](#clustering)
  - [K-means](#k-means)
  - [Hierarchical Clustering](#hierarchical-clustering)
- [Principal Component Analysis
  (PCA)](#principal-component-analysis-pca)
  - [Data import](#data-import)
  - [PCA to the rescue](#pca-to-the-rescue)

Today we will explore unsupervised machine learning methods starting
with clustering and dimensionality reduction.

## Clustering

To start, let’s make up some data to cluster where we know what the
answer should be. The `rnorm()` function will help us here.

``` r
hist(rnorm(10000, mean=3))
```

![](class07_files/figure-commonmark/unnamed-chunk-1-1.png)

Return 30 numbers centered on -3.

``` r
tmp <- c(rnorm(30, mean=-3), rnorm(30, mean=3))

x <- cbind(x=tmp, y=rev(tmp))

x
```

                  x         y
     [1,] -2.256482  3.251890
     [2,] -2.066531  2.309818
     [3,] -2.897704  1.985845
     [4,] -1.343658  1.456434
     [5,] -3.373166  1.262379
     [6,] -2.878685  2.554073
     [7,] -3.063307  4.179239
     [8,] -3.885976  2.957549
     [9,] -4.276082  2.283294
    [10,] -2.699255  2.020119
    [11,] -2.344681  3.696876
    [12,] -3.451357  2.538568
    [13,] -3.239751  2.888730
    [14,] -4.184840  3.629600
    [15,] -2.585861  3.291475
    [16,] -4.214955  1.455182
    [17,] -2.460774  2.650768
    [18,] -1.953049  4.549005
    [19,] -2.618081  2.772930
    [20,] -3.711843  3.335816
    [21,] -2.552060  3.453454
    [22,] -2.640518  1.288053
    [23,] -2.739463  2.065802
    [24,] -2.689124  3.881180
    [25,] -3.655227  2.525409
    [26,] -3.149199  3.731757
    [27,] -2.987928  2.392178
    [28,] -2.618860  5.332014
    [29,] -2.748709  2.652052
    [30,] -2.889492  3.090653
    [31,]  3.090653 -2.889492
    [32,]  2.652052 -2.748709
    [33,]  5.332014 -2.618860
    [34,]  2.392178 -2.987928
    [35,]  3.731757 -3.149199
    [36,]  2.525409 -3.655227
    [37,]  3.881180 -2.689124
    [38,]  2.065802 -2.739463
    [39,]  1.288053 -2.640518
    [40,]  3.453454 -2.552060
    [41,]  3.335816 -3.711843
    [42,]  2.772930 -2.618081
    [43,]  4.549005 -1.953049
    [44,]  2.650768 -2.460774
    [45,]  1.455182 -4.214955
    [46,]  3.291475 -2.585861
    [47,]  3.629600 -4.184840
    [48,]  2.888730 -3.239751
    [49,]  2.538568 -3.451357
    [50,]  3.696876 -2.344681
    [51,]  2.020119 -2.699255
    [52,]  2.283294 -4.276082
    [53,]  2.957549 -3.885976
    [54,]  4.179239 -3.063307
    [55,]  2.554073 -2.878685
    [56,]  1.262379 -3.373166
    [57,]  1.456434 -1.343658
    [58,]  1.985845 -2.897704
    [59,]  2.309818 -2.066531
    [60,]  3.251890 -2.256482

Make a plot of `x`.

``` r
plot(x)
```

![](class07_files/figure-commonmark/unnamed-chunk-3-1.png)

### K-means

The main function in “base” R for K-means clustering is called
`kmeans()`:

``` r
km <- kmeans(x, centers=2)
km
```

    K-means clustering with 2 clusters of sizes 30, 30

    Cluster means:
              x         y
    1 -2.939221  2.849405
    2  2.849405 -2.939221

    Clustering vector:
     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

    Within cluster sum of squares by cluster:
    [1] 40.78012 40.78012
     (between_SS / total_SS =  92.5 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

The `kmeans()` function returns a “list” with 9 components. You can see
the named components of any list with the `attributes()` function.

``` r
attributes(km)
```

    $names
    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

    $class
    [1] "kmeans"

> Q. How many points are in each cluster?

``` r
km$size
```

    [1] 30 30

> Cluster assignment/membership vector?

``` r
km$cluster
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

> Q. Cluster centers?

``` r
km$centers
```

              x         y
    1 -2.939221  2.849405
    2  2.849405 -2.939221

> Q. Make a plot of our `kmeans()` results showing cluster assignment
> using different colors for each cluster/group of points and cluster
> centers in blue.

``` r
plot(x, col=km$cluster)
points(km$centers, col="blue", pch=15, cex=2)
```

![](class07_files/figure-commonmark/unnamed-chunk-9-1.png)

> Q. Run `kmeans()` again on `x` and this cluster into 4 groups/clusters
> and plot the same result figure as above.

``` r
new_km <- kmeans(x, centers=4)
plot(x, col=new_km$cluster)
points(new_km$centers, col="blue", pch=15, cex=2)
```

![](class07_files/figure-commonmark/unnamed-chunk-10-1.png)

> **key-point**: K-means clustering is super popular but can be misused.
> One big limitation is that it can impose a clustering pattern on your
> data even if clear natural grouping don’t exist - i.e. it does what
> you tell it to do in terms of `centers`.

### Hierarchical Clustering

The main function in “base” R for Hierarchical Clustering is called
`hclust()`.

You can’t just pass our dataset as is into `hclust()`. You must give
“distance matrix” as input. We can get this from the `dist()` function
in R.

``` r
d <- dist(x)
hc <- hclust(d)
hc
```


    Call:
    hclust(d = d)

    Cluster method   : complete 
    Distance         : euclidean 
    Number of objects: 60 

The results of `hclust()` don’t have a useful `print()` method, but do
have a special `plot()` method.

``` r
plot(hc)
abline(h=8, col="red")
```

![](class07_files/figure-commonmark/unnamed-chunk-12-1.png)

To get our main cluster assignment (membership vector) we need to “cut”
the tree at the big goal posts…

``` r
grps <- cutree(hc, h=8)
grps
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
table(grps)
```

    grps
     1  2 
    30 30 

``` r
plot(x, col=grps)
```

![](class07_files/figure-commonmark/unnamed-chunk-15-1.png)

Hierarchical Clustering is distinct in that the dendrogram (tree figure)
can reveal the potential grouping in your data (unlike K-means).

## Principal Component Analysis (PCA)

PCA is a common and highly useful dimensionality reduction technique
used in many fields - particularly bioinformatics.

Here we will analyze some data from the UK on food consumption.

### Data import

``` r
url <- "https://tinyurl.com/UK-foods"
x <- read.csv(url)

head(x)
```

                   X England Wales Scotland N.Ireland
    1         Cheese     105   103      103        66
    2  Carcass_meat      245   227      242       267
    3    Other_meat      685   803      750       586
    4           Fish     147   160      122        93
    5 Fats_and_oils      193   235      184       209
    6         Sugars     156   175      147       139

``` r
x <- read.csv(url, row.names=1)
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

``` r
barplot(as.matrix(x), beside=T, col=rainbow(nrow(x)))
```

![](class07_files/figure-commonmark/unnamed-chunk-18-1.png)

``` r
barplot(as.matrix(x), beside=F, col=rainbow(nrow(x)))
```

![](class07_files/figure-commonmark/unnamed-chunk-19-1.png)

One conventional plot that can be useful is called a “pairs” plot.

``` r
pairs(x, col=rainbow(nrow(x)), pch=16)
```

![](class07_files/figure-commonmark/unnamed-chunk-20-1.png)

### PCA to the rescue

The main function in base R for PCA is called `prcomp()`.

``` r
pca <- prcomp(t(x))
summary(pca)
```

    Importance of components:
                                PC1      PC2      PC3       PC4
    Standard deviation     324.1502 212.7478 73.87622 3.176e-14
    Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

The `prcomp()` function returns a list object of our results with five
attributes/components.

``` r
attributes(pca)
```

    $names
    [1] "sdev"     "rotation" "center"   "scale"    "x"       

    $class
    [1] "prcomp"

The two main “results” in here are `pca$x` and `pca$rotation`. The first
of these (`pca$x`) contains the scores of the data on the new PC axis -
we use these to make our “PCA plot”.

``` r
pca$x
```

                     PC1         PC2        PC3           PC4
    England   -144.99315   -2.532999 105.768945 -4.894696e-14
    Wales     -240.52915 -224.646925 -56.475555  5.700024e-13
    Scotland   -91.86934  286.081786 -44.415495 -7.460785e-13
    N.Ireland  477.39164  -58.901862  -4.877895  2.321303e-13

``` r
library(ggplot2)

# Make a plot of pca$x with PC1 vs PC2
# ggrepel , geom_text_repel
ggplot(pca$x) +
  aes(PC1, PC2, label=rownames(pca$x)) +
  geom_point() +
  geom_label()
```

![](class07_files/figure-commonmark/unnamed-chunk-24-1.png)

The graph above shows that the results from N. Ireland are most
different compared to the other three countries.

The second major result is contained in the `pca$rotation` object or
component. Let’s plot this to see what PCA is picking up…

``` r
ggplot(pca$rotation) +
  aes(PC1, rownames(pca$rotation)) +
  geom_col()
```

![](class07_files/figure-commonmark/unnamed-chunk-25-1.png)

The graph above shows that, based on the previous graph, soft drinks and
fresh potatoes are more likely to be associated with N. Ireland and
every other food are more likely to be associated with the other three
countries.
