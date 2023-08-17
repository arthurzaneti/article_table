# For installing
Type the following in your R terminal
```r
devtools::install_github("https://github.com/arthurzaneti/chen")
```

If you get the following message :
```r
Error in loadNamespace(x) : there is no package called ‘devtools’
```
Try dowloading devtools and then dowloading from github, you can copy and paste the following in your R terminal, or Rstudio script:
```r
install.packages("devtools")
devtools::install_github("arthurzaneti/chen")
```
