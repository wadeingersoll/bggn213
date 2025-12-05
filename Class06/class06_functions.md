# Class 6: R Function
Wade Ingersoll (PID: A69038080)

- [Our first (silly) function](#our-first-silly-function)
- [A second function](#a-second-function)
- [Protein generating function](#protein-generating-function)

All functions in R have at least 3 things:

- A **name**, we pick this and use it to call our function,
- Input **arguments** (there can be multiple)
- The **body** lines of R code that do the work

## Our first (silly) function

Write a function to add some numbers

``` r
add <- function(x, y=1) {
  x + y
}
```

Now we can call this function:

``` r
add(c(10, 10), 100)
```

    [1] 110 110

``` r
add(10, 100)
```

    [1] 110

## A second function

Write a function to generate random nucleotide sequences of a user
specified length:

The `sample()` function can be helpful here.

``` r
sample(c("A", "C", "G", "T"), size=50, replace=TRUE)
```

     [1] "A" "G" "A" "T" "G" "A" "C" "C" "G" "G" "G" "A" "G" "G" "G" "C" "A" "G" "T"
    [20] "T" "C" "A" "T" "A" "T" "T" "A" "G" "G" "G" "T" "T" "G" "C" "A" "T" "T" "G"
    [39] "A" "C" "T" "A" "G" "A" "A" "C" "T" "G" "A" "A"

I want a 1-element long character vector that looks like: “CCAGTC”, not
“C” “C” “A” “G” “T” “C”

``` r
v <- sample(c("A", "C", "G", "T"), size=50, replace=TRUE)
paste(v, collapse="")
```

    [1] "CGAGGGTCTTAGGCTGTGTTGAATACCCTTATCACCAATGAATAAGCCCT"

Turn this into my first function

``` r
generate_dna <- function(size = 50) {
  v <- sample(c("A", "C", "G", "T"), size = size, replace = TRUE)
  paste(v, collapse = "")
}
```

Test it:

``` r
generate_dna(30)
```

    [1] "AAGACTTAGGAGTCCCCGCACTGACATGAA"

Practicing IF:

``` r
fasta <- TRUE
if(fasta) {
  cat("HELLO You!")
} else {
  cat("No you don't!")
}
```

    HELLO You!

Add the ability to return a multi-element vector or a single element
fasta-like element:

``` r
generate_dna <- function(size = 50, fasta=TRUE) {
  v <- sample(c("A", "C", "G", "T"), size = size, replace = TRUE)
  s <- paste(v, collapse = "")
  
  if(fasta){
    return(s)
  } else {
    return(v)
  }
}
```

Try it:

``` r
generate_dna(fasta = TRUE)
```

    [1] "TGATGGAACCTTGAGTCTCCTAGCCCGTACCAGCCCCGTGTCCATACGTA"

``` r
generate_dna(fasta = FALSE)
```

     [1] "G" "G" "C" "C" "T" "C" "C" "T" "G" "A" "C" "C" "C" "T" "T" "T" "A" "T" "T"
    [20] "A" "T" "C" "G" "C" "C" "A" "T" "C" "G" "A" "G" "A" "T" "A" "C" "C" "T" "A"
    [39] "C" "A" "G" "C" "G" "A" "G" "A" "C" "T" "T" "A"

## Protein generating function

``` r
generate_protein <- function(size = 50, fasta = TRUE) {
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", 
          "H", "I", "L", "K", "M", "F", "P", "S", 
          "T", "W", "Y", "V")
  
  v <- sample(aa, size = size, replace = TRUE)
  s <- paste(v, collapse = "")
  
  if (fasta) {
    return(s)
  } else {
    return(v)
  }
}
```

Try it:

``` r
generate_protein(10)
```

    [1] "HTKMRMAEIK"

Use our new `generate_protein()` function to make random protein
sequences of length 6 to 12 (i.e. one length 6, one length 7, etc. up to
length 12)

One way to do this is “brute force”

``` r
generate_protein(6)
```

    [1] "SYTDKD"

``` r
generate_protein(7)
```

    [1] "SCGCPFQ"

``` r
generate_protein(8)
```

    [1] "ANFRNGQP"

``` r
generate_protein(9)
```

    [1] "HHLMKIWWC"

A second way is to use a `for()` loop:

``` r
lengths <- 6:12
lengths
```

    [1]  6  7  8  9 10 11 12

``` r
for(i in lengths) {
  cat(">", i, '\n', sep="")
  aa <- generate_protein(i)
  cat(aa)
  cat('\n')
}
```

    >6
    VIQPTW
    >7
    SMGCRAT
    >8
    NITPPYHQ
    >9
    DITIYDLQF
    >10
    INTHEYSNHG
    >11
    ICQICWRRYMC
    >12
    YFMRFMSGRKGG

Try `sapply()` function

``` r
sapply(6:12, generate_protein)
```

    [1] "HAEPRD"       "DSQCDYF"      "SMQVCTFW"     "VGPPNQAFA"    "NYMTEEKYDY"  
    [6] "ATRYCGTIDWF"  "WFMMGFQGCAIY"
