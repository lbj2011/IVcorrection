# IV curve correction tools (IVcorrection)

<img src="https://github.com/lbj2011/IVcorrection/blob/main/doc_img/ivcorrection_logo.png" width="500"/>


[![PyPI release](https://img.shields.io/pypi/v/ivcorrection.svg)](https://pypi.org/project/ivcorrection/)
[![DOI](https://zenodo.org/badge/681861201.svg)](https://zenodo.org/badge/latestdoi/681861201)

This package provides the following purposes:
 - Calculate IV curve correction coefficients based on IEC 60891:2021 standard
 - Perform IV curve correction

**An interactive page for this tool will be online soon!**

**The package is still under active development and plan to be integrated with [pvlib](https://github.com/pvlib/pvlib-python).**
 
**If there is any problem, please feel free to [contact us](mailto:baojieli@lbl.gov)!**

## Installation

- Install with pip
```
pip install ivcorrection
```

## What is I-V curve correction
I-V curve correction refers to correct the I-V curve measured under different environmental condition to an identical one. Generally, the standard test condition (STC) is adopted as the target condition.

## Correction methods
This tool can perform I-V curve correction using **Procedure 1, 2, 4 of IEC 60891:2021 [1]** and a **proposed Procedure dynamic**. The determination of correction coefficients follow the IEC 60891:2021 standard [1].

## How to use this tool
The calculation proceeds with the following steps:

- Provide PV module parameters
- Calculate correction coefficients
- Correct example I-V curve(s) (optional)

An example is provided in [Example_IV_correction_CEC.ipynb](https://github.com/lbj2011/IVcorrection/tree/main/examples/Example_IV_correction_CEC.ipynb).
<img src="https://github.com/lbj2011/IVcorrection/blob/main/doc_img/example.png" width="700"/>


## Who we are
We are a collection of national lab researchers funded under the Durable module materials consortium (DuraMAT).
<img src="https://github.com/lbj2011/IVcorrection/blob/main/doc_img/duramat_logo.png" width="200"/>

## References
[1] IEC 60891, Photovoltaic devices - Procedures for temperature and irradiance corrections to measured I-V characteristics, 2021.


## Authors
Baojie Li (LBNL)