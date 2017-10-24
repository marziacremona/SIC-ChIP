# SIC-ChIP
#### Shape Index Clustering for ChIP-seq peaks

<img src="https://raw.githubusercontent.com/marziacremona/SIC-ChIP/master/SIC-ChIP_logo.png" alt="SIC-ChIP logo" width="200">

In the last years, ChIP-seq (Chromatin ImmunoPrecipitation sequencing) experiments have been widely used to investigate protein-DNA interactions. At present many algorithms and pipelines are available for downstream analysis of these datasets, each looking in a slightly different way at the intensity of ChIP-seq signals to detect the enriched regions (peaks). Although peak shape can vary widely, both among ChIP-seqs for different proteins and among different enriched regions in a single ChIP-seq, it is almost always ignored by these analysis techniques. 

We propose a novel analysis pipeline that takes into consideration both the shape and the intensity of the peaks in a single ChIP-seq, with the idea that peak shape might provide scientists with new insights into chromatin function. Specifically, we select five indices of shapes (that are particularly suited for transcription factor peaks) and we use them to cluster ChIP-seq peaks by looking at their shapes. 

SIC-ChIP (Shape Indices Clustering for ChIP-seq peaks) implements the core of our methodology. Given the ChIP-seq coverage function (BigWig file) and the list of peaks (BED file), SIC-ChIP defines the peak curves, compute the five indices of shape and uses k-mean algorithm to cluster the peaks, testing different numbers of clusters. 

### Publication

Detailed information about SIC-ChIP and the complete analysis pipeline can be found in [Cremona, Sangalli, Vantini, Dellino, Pelicci, Secchi, Riva. **Peak shape clustering reveals biological insights**. 2015. **16**:349.](https://doi.org/10.1186/s12859-015-0787-6)


