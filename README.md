# Kubota_puQTL

This repository contains the data analysis code used in the manuscript:

**Mapping of promoter usage QTL using RNA-seq data reveals their contributions to complex traits**
[https://www.biorxiv.org/content/10.1101/2022.02.24.481875v1](https://www.biorxiv.org/content/10.1101/2022.02.24.481875v1)

## Dependencies

- Python (3.8.5)
- Pandas (1.1.3)
- Matplotlib (3.3.1)
- Seaborn (0.11.0)
- Scipy (1.6.2)
- Docker (20.10.7)

### Docker images
- [naotokubota/proactiv:1.1.18](https://hub.docker.com/repository/docker/naotokubota/proactiv)
- [broadinstitute/gtex_eqtl:V8](https://hub.docker.com/r/broadinstitute/gtex_eqtl)
- [naotokubota/qtltools:1.3.1](https://hub.docker.com/repository/docker/naotokubota/qtltools)
- [naotokubota/coloc-locuscomparer:1.0](https://hub.docker.com/repository/docker/naotokubota/coloc-locuscomparer)
- [nservant/hicpro:3.0.0](https://hub.docker.com/r/nservant/hicpro)
- [aylab/fithichip:latest](https://hub.docker.com/r/aylab/fithichip)
- [quay.io/biocontainers/deeptools:3.5.1--py_0](https://quay.io/repository/biocontainers/deeptools)
- [quay.io/biocontainers/plink:1.90b6.21--h779adbc_1](https://quay.io/repository/biocontainers/plink)