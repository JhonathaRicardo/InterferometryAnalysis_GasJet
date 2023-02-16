# Interferometry Analysis - Gas-Jet Profile

## Abstract
The interferometric technique is an important tool for analysis and diagnosis in astronomy, spectroscopy, metrology, plasma physics, particle physics, and other areas. In Laser Wakefield Acceleration (LWFA) studies, knowing the density distribution of the target gas is crucial to understand the phenomena involved in the particle acceleration process [[1]](#reference).
This Python algorithm was developed to recover the accumulated phase across the neutral gas target as well as estimate the target density distribution.

#### Badges
![License](https://img.shields.io/badge/license-MIT-green)
![version](https://img.shields.io/badge/version-v.1.0-green)
![status](https://img.shields.io/badge/status-under%20development-green)
## Summary
* [Abstract](#abstract)
* [Badges](#badges)
* [Summary](#summary)
* [Installation](#installation)
* [Authors](#pessoas-desenvolvedoras)
* [License](#license)
* [Citation](#citation)
## Installation

The Interferometry Analysis - Gas-Jet software was developed in Python 3.11 and the use of this algorithm requires the installation of some packages: [NumPy](https://numpy.org/), [Scipy](https://scipy.org/) and [PyAbel](https://pyabel.readthedocs.io/en/latest/index.html) for data processing, [Pillow](https://pypi.org/project/Pillow/) to manipulate interferogram images, [Matplotlib](https://matplotlib.org/stable/index.html) to plot results, and
[PySimpleGui](https://www.pysimplegui.org/en/latest/) to create the users template.

The second way to use this software is through the executable file. The users can create a single .exe file using the [pyinstaller](https://pyinstaller.org/en/stable/) package trought the follow terminal command:

<code>   pyinstaller --onefile -w intAnalysis_GasJetProfile.py                </code>

## Reference
[1] Hariharan, P. (2007) Basics of Interferometry. 2nd Edition, Elsevier, Amsterdam.[https://doi.org/10.1016/B978-0-12-373589-8.X5000-7](https://doi.org/10.1016/B978-0-12-373589-8.X5000-7) 
## License
## Citation

