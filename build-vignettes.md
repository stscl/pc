Due to the time-consuming computations involved in the vignettes of the *pc* package, 
it is necessary to pre-build the vignettes prior to package submission.

``` r
.prebuild_vignettes = \(name){
  out = paste0("vignettes/",name,".Rmd")
  inp = paste0(out,".orig")
  knitr::knit(inp,out)
}

.prebuild_vignettes("pc")
```
