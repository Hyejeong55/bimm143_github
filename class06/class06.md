# Class 6: R Functions
Hyejeong Choi (PID: A16837133)

- [1. Generating functions](#1-generating-functions)
- [2. Generate DNA sequence](#2-generate-dna-sequence)
- [3. Generate protein sequence](#3-generate-protein-sequence)

## 1. Generating functions

Let’s start writing our first silly function to add some numbers:

Every R function has 3 things:

- name (we get to pick this)
- input arguments (there can be loads of these separated by a comma)
- the body (the R code that does the work)

``` r
add <- function(x, y=10, z=0){x + y + z}
```

I can just use this function like any other function as long as R knows
about it (i.e. run the code chunk)

``` r
add(1, 100)
```

    [1] 101

``` r
add(x=c(1,2,3,4), y=100)
```

    [1] 101 102 103 104

``` r
add(1)
```

    [1] 11

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with an equals default
value (`y=0`) in the function definition.

``` r
add(1, 100, 10)
```

    [1] 111

> Q. Write a function to return a DNA sequence of a user specified
> length? Call it `generate_dna()`

The `sample()` function can help here

``` r
# generate_dna <- function(size=5){}

students <- c("jeff", "jeremy", "peter")

sample(students, size = 5, replace = TRUE)
```

    [1] "jeremy" "peter"  "jeremy" "jeff"   "peter" 

## 2. Generate DNA sequence

Now work with `bases` rather than `students`

``` r
bases <- c("A", "C", "G", "T")

sample(bases, size = 10, replace = TRUE)
```

     [1] "G" "A" "G" "G" "C" "A" "G" "T" "T" "T"

Now I have a working ‘snippet’ of code I can use as the body of my first
function version here:

``` r
generate_dna <- function(length=5){
  bases <- c("A", "C", "G", "T")
  sample(bases, size=length, replace=TRUE)}
```

``` r
generate_dna()
```

    [1] "A" "T" "A" "A" "C"

I want the ability to return a sequence like “AGTACCTG” i.e. a one
element vector where the bases are all together.

``` r
generate_dna <- function(length=5, together=TRUE){
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, size=length, replace = TRUE)
  
  if(together){
    sequence <- paste(sequence, collapse = "")}
  
  return(sequence)}
```

``` r
generate_dna()
```

    [1] "CGACA"

## 3. Generate protein sequence

> Q. Write a protein sequence generating function that will return
> sequences of a user specified length?

``` r
generate_protein <- function(length=10, together=TRUE){
  
  ## get the 20 amino acids as a vector
  amino_acids <- bio3d::aa.table$aa1[1:20]
  sequence_protein <- sample(amino_acids, size=length, replace=TRUE)
  
  ## optionally return a single element string
  if(together){
    sequence_protein <- paste(sequence_protein, collapse="")}
  
  return(sequence_protein)}
```

``` r
generate_protein()
```

    [1] "SRIVVNCGQY"

> Q. Generate random protein sequences of length 6 to 12 amino acids.

``` r
generate_protein()
```

    [1] "QGFITRGNKQ"

We can fix this inability to generate multiple sequences by either
editing and adding to the function body code (e.g. a for loop) or by
using the R **apply** family of utility functions.

``` r
sapply(6:12, generate_protein)
```

    [1] "WMDNTY"       "PSRWICY"      "WYYPLRQW"     "IHSPMQFYN"    "WMFAMQITNI"  
    [6] "RDYDPGTRIFL"  "DWTMFTYNTIHR"

It would be cool and useful if I could get FASTA format output.

``` r
ans <- sapply(6:12, generate_protein)
ans
```

    [1] "CCESAG"       "YPTPCGT"      "GKSEADPR"     "LEQTEMTCM"    "ADCFTKFVCA"  
    [6] "YMTKCWGVQFF"  "LTSRLYWGGYQN"

``` r
cat(ans, sep="\n")
```

    CCESAG
    YPTPCGT
    GKSEADPR
    LEQTEMTCM
    ADCFTKFVCA
    YMTKCWGVQFF
    LTSRLYWGGYQN

I want this to look like:

> ID.6 KIWHND ID.7 FENRCVK ID.8 MSKHTDFI

The functions `paste()` and `cat()` can help us here

``` r
cat(paste("> ID.", 6:12, "\n", ans, sep=""), sep = "\n" )
```

    > ID.6
    CCESAG
    > ID.7
    YPTPCGT
    > ID.8
    GKSEADPR
    > ID.9
    LEQTEMTCM
    > ID.10
    ADCFTKFVCA
    > ID.11
    YMTKCWGVQFF
    > ID.12
    LTSRLYWGGYQN

``` r
id.line <- paste("> ID.", 6:12, sep="")
id.line
```

    [1] "> ID.6"  "> ID.7"  "> ID.8"  "> ID.9"  "> ID.10" "> ID.11" "> ID.12"

``` r
seq.line <- paste(id.line, ans, sep="\n")
cat(seq.line, sep="\n", file="myseq.fa")
```

> Q. Determine if these sequences can be found in nature or are they
> unique? Why or why not?

I BLASTp searched my FASTA format sequences against the refseq_protein
database and found that length 6, 7, and 8 are not unique and can be
found in the database with 100% coverage and 100% identity. Length 9,
10, 11, and 12 are unique and can’t be found in the database.

We can get the set of 20 natural amino acids from the **bio3d** package.

``` r
bio3d::aa.table$aa1[1:20]
```

     [1] "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T" "W" "Y"
    [20] "V"
