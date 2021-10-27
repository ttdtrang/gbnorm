# gbnorm - A graph-based algorithm for RNA-seq data normalization

Citation:
* Tran DT, Bhaskara A, Kuberan B, Might M (2020) A graph-based algorithm for RNA-seq data normalization. PLOS ONE 15(1): e0227760. https://doi.org/10.1371/journal.pone.0227760


## Notes

* [Release v0.99.1](https://github.com/ttdtrang/gbnorm/releases/tag/v0.99.1) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3628859.svg)](https://doi.org/10.5281/zenodo.3628859) is the first working snapshot of gbnorm. This version was used in the paper https://doi.org/10.1371/journal.pone.0227760. In this implementation, all the functions taking expression matrix as an input assume that it is _samples_ (rows) x _features_ (columns).
* [Release v0.99.2](https://github.com/ttdtrang/gbnorm/releases/tag/v0.99.2) made changes to the major functions that take expression matrix as an input such that the matrix is _features_ (rows) x _samples_ (columns).

## License

GNU GPLv3
