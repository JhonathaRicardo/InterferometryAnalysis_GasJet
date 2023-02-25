# <h1 align = "center">Interferometry Analysis - Gas-Jet Profile </h1>
![Logo](https://github.com/JhonathaRicardo/InterferometryAnalysis_GasJetProfile/blob/main/Images/Intro_GasJet.jpg)
![License](https://img.shields.io/badge/license-MIT-green)
![version](https://img.shields.io/badge/version-v.1.0-green)
![status](https://img.shields.io/badge/status-under%20development-green)
## Abstract
The interferometric technique is an important tool for analysis and diagnosis in astronomy, spectroscopy, metrology, plasma physics, particle physics, and other areas. In Laser Wakefield Acceleration (LWFA) studies, knowing the density distribution of the target gas is crucial to understand the phenomena involved in the particle acceleration process [[1]](#reference).
This Python algorithm was developed to recover the accumulated phase across the neutral gas target as well as estimate the target density distribution.

## Summary
* [Installation](#installation)
* [How to use](#how-to-use)
  * [Main Screen](#main-screen)
* [Reference](#reference)
* [Authors](#authors)
* [License](#license)
* [Citation](#citation)

## Installation
Interferometry Analysis - Gas-Jet software was developed in Python 3.11 and the use of this algorithm requires the installation of some packages: [NumPy](https://numpy.org/), [Scipy](https://scipy.org/) and [PyAbel](https://pyabel.readthedocs.io/en/latest/index.html) for data processing, [Pillow](https://pypi.org/project/Pillow/) to manipulate interferogram images, [Matplotlib](https://matplotlib.org/stable/index.html) to plot results, and
[PySimpleGui](https://www.pysimplegui.org/en/latest/) to create the users template.

The second way to use this software is through the executable file. The users can create a single .exe file using the [pyinstaller](https://pyinstaller.org/en/stable/) package trought the follow terminal command:

<code>   pyinstaller --onefile -w intAnalysis_GasJetProfile.py                </code>

## How to use
Interferometry Analysis - The Gas-Jet software has a graphical interface developed with PysimpleGUI. This interface assists users and facilitates their applications.
In this section, we provide users with a simple review of the software's functions and how to use them.
### Main Screen
![Template](https://github.com/JhonathaRicardo/InterferometryAnalysis_GasJetProfile/blob/main/Images/Template_image.jpg)
### Interferogram Image Frame
- ***[Interferogram (Gas-Jet)]*** Scaled gas-jet interferogram image.

- ***[Open File(s)]*** Open interferogram image(s) file(s) with the presence of gas jet . Image file extensions should preferably be .png or .snp. However, all image extensions (.gif, .jpg, .bmp, etc) could be used. The path to opened file is shown in text box above.  

- ***[Interferogram (Ref.)]*** Scaled reference interferogram image.

- ***[Open File]*** Open an undisturbed interferogram image(s) file(s). Image file extensions should preferably be .png or .snp. However, all image extensions (.gif, .jpg, .bmp, etc) could be used. The path to opened file is shown in text box above.  
 
- ***[Rotate]*** The image rotates in degrees. Positive degrees promote counterclockwise rotation.  

- ***[Image Scale]*** The interferogram image shown is scaled to screen size (428,342) for users' viewing only. However, all processes to determine the gas jet density profile are done with the original dimensions of the image file.

- ***[Analyse Data]*** From this command button, the software will apply data processing to generate accumulated phase, inverse Abel transforms, and gas jet density profile.

- ***[Exit]*** Close software.

### Options
- ***[Select Analysis Area]*** From the parameters in this form, the user can select the interferogram area to apply the algorithm to determine the gas-jet density profile. The selected area is defined by a rectangle with edges defined by X and Y coordinates(***[Y Coord]*** and ***[X Coord]***).
The user that intends to use the whole interferogram figure needs to uncheck the checkbox ***[Use select area]***.

- ***[Experimental Parameteres]*** In this form, the user set the experimental parameters used to obtain the interferogram image. These parameters are:
  - ***[Laser Wavelength]*** and ***[uncertainty Laser Wavelength]*** in nm;
  - Gas ***[Polarizability]*** in angstrom³. This parameter usually refers to the tendency of matter to acquire an electric dipole moment when subjected to an electric field.
   - The ***[Specific Heat Ratio]*** of a gas is the ratio of the specific heat at constant pressure, ***Cp***, to the specific heat at constant volume, ***Cv***.

- ***[Analysis Parameters]*** This form contains the parameters for analysis of the interferogram images:
  - ***[Scaling Factor]*** of an interferogram image in pixels/micrometers;
  - ***[Sigma - Gaussian Blur]*** is a multidimensional gaussian image filter. The standard deviation of the gaussian filter (Sigma) defined by the user is equal for all axes.
  - ***[ Gaussian Filter Position]*** this parameter is different from Gaussian Blur. This parameter is set automatically by the algorithm and this position defines which frequency will be used to apply the Inverse Fourier Transform and build the phase map of the gas-jet.
Both the above parameters are defined in pixels. ***Note:*** The algorithm set the frequency that defines a positive phase map.  But, users can change the filter position.
  - ***[Fringes Orientation]*** can be vertical or horizontal.
  - ***[Axisymmetric]*** An important parameter to apply the Inverse Abel Transform is the axis of symmetry (or axisymmetric). The axisymmetric can be horizontal or vertical.

### Gas-Jet Profile
- ***[Stages]:*** The stages of the results obtained by the algorithm can be viewed by user.
 - ***[Fourier Transform]*** This image is build trought the Fourier Transform of interferogram image.


## Reference
[1] Hariharan, P. (2007) Basics of Interferometry. 2nd Edition, Elsevier, Amsterdam.[https://doi.org/10.1016/B978-0-12-373589-8.X5000-7](https://doi.org/10.1016/B978-0-12-373589-8.X5000-7)
## Authors
Interferometry Analysis - Gas-Jet software was developed by researchs of the High Power Ultrashort Pulse Lasers Group of the Center for Lasers and Applications at [IPEN](https://www.ipen.br/portal_por/portal/default.php) (Instituto de Pesquisas Energéticas e Nucleares).
* Jhonatha Ricardo dos Santos [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7877-0580)
* Armando Valter Felicio Zuffi [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-5705-1499)
* Ricardo Elgul Samad [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7762-8961)
* Nilson Dias Vieira Junior [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0003-0092-9357)

## License
Copyright (c) 2023 Jhonatha Ricardo dos Santos

## Citation

