# SIC-ChIP - Shape Index Clustering for ChIP-seq peaks

<img src="https://raw.githubusercontent.com/marziacremona/SIC-ChIP/master/SIC-ChIP_logo.png" alt="SIC-ChIP logo" width="200">

In the last years, ChIP-seq (Chromatin ImmunoPrecipitation sequencing) experiments have been widely used to investigate protein-DNA interactions. At present many algorithms and pipelines are available for downstream analysis of these datasets, each looking in a slightly different way at the intensity of ChIP-seq signals to detect the enriched regions (peaks). Although peak shape can vary widely, both among ChIP-seqs for different proteins and among different enriched regions in a single ChIP-seq, it is almost always ignored by these analysis techniques. 

We propose a novel analysis pipeline that takes into consideration both the shape and the intensity of the peaks in a single ChIP-seq, with the idea that peak shape might provide scientists with new insights into chromatin function. Specifically, we select five indices of shapes (that are particularly suited for transcription factor peaks) and we use them to cluster ChIP-seq peaks by looking at their shapes. 

SIC-ChIP (Shape Indices Clustering for ChIP-seq peaks) implements the core of our methodology. Given the ChIP-seq coverage function (BigWig file) and the list of peaks (BED file), SIC-ChIP defines the peak curves, compute the five indices of shape and uses k-mean algorithm to cluster the peaks, testing different numbers of clusters. 

## Publication
Detailed information about SIC-ChIP and the complete analysis pipeline can be found in [Cremona, Sangalli, Vantini, Dellino, Pelicci, Secchi, Riva. **Peak shape clustering reveals biological insights**. 2015. **16**:349](https://doi.org/10.1186/s12859-015-0787-6).

## Documentation
SIC-ChIP works on UNIX and MAC systems.

### Requirements: 
- R (version >3.1.2)
- R packages [*rtracklayer*](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html), [*parallel*](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/parallel-package.html), [*BH*](https://cran.r-project.org/web/packages/BH/index.html) and [*Rcpp*](https://cran.r-project.org/web/packages/Rcpp/index.html)
- [*Boost*](http://www.boost.org) C++ libraries

### Usage:
- Open a new terminal and change directory to the folder where your input files *coverage_test.bw* and *peaks_test.bed* are located
- Run SIC-ChIP with the following command

  ```/path/to/SIC-ChIP.r --bw=coverage_test.bw --bed=peaks_test.bed --out=results [options]```

  or alternatively: 
  
  ```Rscript /path/to/SIC-ChIP.r --bw=coverage_test.bw --bed=peaks_test.bed --out=results [options]```

##### Mandatory options:
   ```--bw``` BigWig file with the coverage function
   
   ```--bed``` BED file with the coordinates of the peaks
   
   ```--out``` The name of the results folder/files

##### Options:
   ```--help``` print short help message and exit
   
   ```--N``` the distance between spline knots for computing the number of local maxima (default ```--N=20```) 
   
   ```--toll``` the minimum distance allowed between local maxima (default ```--toll=50```)

**IMPORTANT**: the parameters ```N``` and ```toll``` should be carefully chosen depending on the peaks dataset (for example, they should be lower than the default when some peaks have low width)


### Output:
A folder named SIC-ChIP_results containing:

- The RData file with a *RangedData* of all the peaks locations and counts, i.e. the coverage function for each peak (*results_peaks.RData*)
- The RData file with a dataframe of the five shape indices for all the peaks (*results_peaks_index.RData*)
- The RData file with a list *K* of k-means results for clusters number *k=1,...10* and a vector *within* with the total within clusters sum of squares plot for cluster number *k=1,...,10* (*results_kmeans.RData*)
- The scatterplot of the five shape indices for all the peaks (*results_scatterplot.pdf*)
- The total within clusters sum of squares plot (*results_tot_within_ss.pdf*)

