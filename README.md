[![test](https://github.com/wtsi-hpag/LinkStats/actions/workflows/test.yml/badge.svg)](https://github.com/wtsi-hpag/LinkStats/actions/workflows/test.yml)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/linkstats/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/linkstats/badges/downloads.svg)](https://anaconda.org/bioconda/linkstats)
# LinkStats
Collect and process statistics from aligned linked-reads.

# Bioconda
LinkStats is available on [bioconda](https://bioconda.github.io/).<br/>
```sh
> conda install linkstats
```

# Usage
```bash
LinkStats --help
```

# Requirments, Running
* [htslib](https://www.htslib.org/) >=1.12
* [python](https://www.python.org/) >=3.8
    * [click](https://click.palletsprojects.com/en/8.0.x/) >=8.0.1
    * [pandas](https://pandas.pydata.org/) >=1.2.5
    * [numpy](https://numpy.org/) >=1.20.3
    * [seaborn](https://seaborn.pydata.org/) >=0.11.1
    * [matplotlib](https://matplotlib.org/stable/index.html) >=3.4.2
    * [tqdm](https://tqdm.github.io/) >=4.61.1

# Requirments, Installing
* [meson](https://mesonbuild.com/)
* C++ compiler (tested with [clang](https://clang.llvm.org/) 11)
* [setuptools](https://setuptools.readthedocs.io/en/latest/)
* [git](https://git-scm.com/)

```bash
cd LinkStats
python setup.py install
```

# Third-Party Acknowledgements
* [stb_sprintf.h](https://github.com/nothings/stb/blob/master/stb_sprintf.h)
