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
  * [Interferogram Image](#interferogram-image)
  * [Options](#options)
  * [Gas-Jet Profile](#gas-jet-profile)
* [How to work](#how-to-work)
* [Reference](#reference)
* [Authors](#authors)
* [Acknowledgment](#acknowledgment)
* [License](#license)
* [Citation](#citation)

## Installation
Interferometry Analysis - Gas-Jet software was developed in Python 3.11 and the use of this algorithm requires the installation of some packages: [NumPy](https://numpy.org/) [[2]](#reference), [Scipy](https://scipy.org/) [[3]](#reference) and [PyAbel](https://pyabel.readthedocs.io/en/latest/index.html) [[4]](#reference) for data processing, [Pillow](https://pypi.org/project/Pillow/) [[5]](#reference)to manipulate interferogram images, [Matplotlib](https://matplotlib.org/stable/index.html) [[5]] to plot results, and
[PySimpleGui](https://www.pysimplegui.org/en/latest/) to create the users template.

The second way to use this software is through the executable file. The users can create a single .exe file using the [pyinstaller](https://pyinstaller.org/en/stable/) package trought the follow terminal command:

<code>   pyinstaller --onefile -w IntAnalysis_GasJetProfile.py                </code>

## How to use
Interferometry Analysis - The Gas-Jet software has a graphical interface developed with PysimpleGUI. This interface assists users and facilitates their applications.
In this section, we provide users with a simple review of the software's functions and how to use them.
### Main Screen

| ![Template](/Images/Template_image.jpg) |
|:--:| 
| *Fig.1 - Software Template* |

### Interferogram Image
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
  - ***[Fourier Transform]*** This image is built through the Fourier Transform of gas-jet interferogram image. From this frequency map, the software selects automatically the frequency that generates a positive phase map. The pixel position (red line) of the selected frequency is the ***[Gaussian Filter position]***. 

    | ![Fourier map](/Images/Fourier_transform_map.png) |
    |:--:| 
    | *Fig.2 - Fourier Transform of gas-jet interferogram* |
  
  > **Note**   
  >  Case the ***[Gaussian Filter position]*** is zero, the software will set the valor automatically.  The user can change this ***[Gaussian Filter position]*** manually.
  - ***[Gaussian Filter]*** This image is the Gaussian filter map applied to generate the phase map using the selected frequency.

    | ![Gaussian filter](/Images/Gaussian_Filter.png) |
    |:--:| 
    | *Fig.3 -  Applied Gaussian filter* |

From the next three steps, users have the option of viewing the average values of the maps in 2D (using ***[2D Profile]*** button) or the 1D profile of the maps in different positions on a symmetrical axis (using ***[1D Profile]*** button). 2D standard deviation maps or 1D curves can be viewed using the ***[Standard Deviation]*** selection box.
  
  - ***[Accumulated Phase]*** Accumulated phase map of the gas-jet.
  
    | ![Acc.phase map](/Images/Accumulated_Phase.png) |
    |:--:| 
    | *Fig.4 -  Accumulated phase of gas-jet interferogram* |   
    
  - ***[Abel Transform]*** Phase map obtained after applying Inverse Abel Transform at the Accumulated Phase map.
  
    | ![Phase map](/Images/Abel_Transform.png) |
    |:--:| 
    | *Fig.5 -  Phase map of gas-jet interferogram* |  
    
  - ***[Density Profile]*** Phase map obtained after applying Inverse Abel Transform at the Accumulated Phase map.
    
    | ![Gasjet_density](/Images/GasjetDensity.png) |
    |:--:| 
    | *Fig.6 -  Gas-jet density profile* |
    
    
## How to work
The interferogram analysis software algorithm works according to the flowchart below:

### Accumulated Phase

#### Standard Deviation of Accumulated Phase

### Inverse Abel Transform

#### Standard Deviation of Inverse Abel Transform

### Density Profile

#### Standard Deviation of Density
  
## Reference
- [1] Hariharan, P. (2007) Basics of Interferometry. 2nd Edition, Elsevier, Amsterdam.[DOI: 10.1016/B978-0-12-373589-8.X5000-7](https://doi.org/10.1016/B978-0-12-373589-8.X5000-7).
- [2] Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). [DOI: 10.1038/s41586-020-2649-2](https://www.nature.com/articles/s41586-020-2649-2). 
- [3] Pauli Virtanen, et. al. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272. [DOI: 10.1038/s41592-019-0686-2](https://www.nature.com/articles/s41592-019-0686-2).
- [4] Gibson, Stephen; Hickstein, Daniel D.; Yurchak, Roman; Ryazanov, Mikhail; Das, Dhrubajyoti; Shih, Gilbert.(2022) PyAbel, PyAbel: v0.9.0, Zenodo,  [DOI: 10.5281/zenodo.7438595](https://doi.org/10.5281/zenodo.7438595).
- [5] Clark, A. (2015). Pillow (PIL Fork) Documentation. readthedocs. Retrieved from [https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf](https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf).
- [6] 

## Authors
Interferometry Analysis - Gas-Jet software was developed by researchs of the High Power Ultrashort Pulse Lasers Group of the Center for Lasers and Applications (CLA) at Instituto de Pesquisas Energéticas e Nucleares ([IPEN](https://www.ipen.br/portal_por/portal/default.php)).
* Jhonatha Ricardo dos Santos [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7877-0580)
* Armando Valter Felicio Zuffi [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-5705-1499)
* Ricardo Elgul Samad [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7762-8961)
* Nilson Dias Vieira Junior [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0003-0092-9357)

## Acknowledgment
Interferogram Analysis - Gas-Jet software was developed to help with the analyze the density of supersonic gas-jet generated by micro-nozzles. This software is part of the work realized by the High Power Ultrashort Pulse Lasers Group of the CLA/IPEN in the [Extreme Light Laboratory](https://www.unl.edu/diocles/home) (ELL)  at the University of Nebraska - Lincoln (UNL). This partnership was able due to funding provided by the São Paulo Research Foundation (FAPESP) and by the US Department of Energy, through the LaserNetUS program.

The author Jhonatha Ricardo dos Santos also acknowledges the FAPESP for doctoral fellowship 2017/13737-8. 

## License
Copyright (c) 2023 Jhonatha Ricardo dos Santos

## Citation

