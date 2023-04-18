# <h1 align = "center">Interferometry Analysis - Gas-Jet (v.1.2)</h1>
<p align="justify">
  Interferometric techniques are important tools for analysis and diagnosis in astronomy, spectroscopy, metrology, plasma physics, particle physics, and other areas, frequently applied to quantify changes in the refractive index of a material or a medium. For example, knowing the density distribution of a gas target is crucial to understand laser plasmas interactions and processes. This software was developed in Python to recover the accumulated optical phase across a neutral gas target, as well as estimate the target density distribution.
</p>

<p align="center">
  <img src = '/Images/Intro_GasJet.jpg' width="80%" align="center">
</p>

![License](https://img.shields.io/badge/license-MIT-green)
![version](https://img.shields.io/badge/version-v.1.2-green)
![status](https://img.shields.io/badge/status-under%20development-yellow)
![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.7807099-blue)

## Summary
* [Introduction](#introduction)
* [Installation](#installation)
* [How to use it](#how-to-use-it)
  * [Main Screen](#main-screen)
  * [Interferogram Image](#interferogram-image)
  * [Options](#options)
  * [Gas-Jet Profile](#gas-jet-profile)
* [How it works](#how-it-works)
* [Example](#example)
* [Authors](#authors)
* [Acknowledgment](#acknowledgment)
* [License](#license)
* [Citation](#citation)
* [Reference](#reference)

## Introduction
  The development of diagnostic tools is very important for a better understanding of laser-plasma interactions [[1]](#reference). An accurate diagnostic is crucial, as instabilities in both target and laser pulses can result in low reproducibility of processes and impair the quality of the intended interaction [[2]](#reference). Among the various non-perturbing optical methods that can be used to diagnose a gaseous target [[3-6]](#reference), interferometry is a very accurate technique capable of quantifying very small optical path differences and, therefore, suitable for measuring density variations of gases [[7, 8]](#reference) and laser-induced plasmas [[1]](#reference). The main drawback of the technique is that returns the integrated phase along the light path, requiring deconvolution methods for retrieving the target density profile. The software “Interferometry Analysis – Gas-Jet” was developed due to the need for a new diagnostic tool to aid in the characterization of supersonic gas jets, quickly and reliably. It was developed by our research group as part of the work to implement a laser-plasma accelerator infrastructure at the Nuclear and Energy Research Institute (IPEN), in Brazil.

## Installation
The *Interferometry Analysis - Gas-Jet* software was developed in Python 3.11. The use of this algorithm requires the installation of the following packages: [NumPy](https://numpy.org/) [[9]](#reference), [Scipy](https://scipy.org/) [[10]](#reference) and [PyAbel](https://pyabel.readthedocs.io/en/latest/index.html) [[11]](#reference) for data processing, [Pillow](https://pypi.org/project/Pillow/) [[12]](#reference) and Scikit-image [[13]](#reference) for the procrssing of interferogram images, [Matplotlib](https://matplotlib.org/stable/index.html) [[14]](#reference) to plot results, and [PySimpleGui](https://www.pysimplegui.org/en/latest/) to create the user's template.

Users also can create a single .exe file using the [pyinstaller](https://pyinstaller.org/en/stable/) package trought the follow terminal command:

<code>   pyinstaller --onefile -w IntAnalysis_GasJet_v1.2.py                </code>

Users who do not use Python IDEs can utilize the software through the executable file available for download [here](https://drive.google.com/file/d/1K_-3wm8TOzxROLyAc22AkRXlRvD3GjjI/view?usp=sharing)

## How to use it
The “Interferometry Analysis – Gas-Jet” has a graphical interface to facilitate its use, and this section provide a simple review of the software's functions and how to employ them.

### Main Screen

|<img src = '/Images/MainScreen1.png'> |
|:--:| 
| *Fig.1. Software Main Screen* |

### Interferogram Image
- ***[Interferogram (Gas-Jet)]*** interferogram image frame.

- ***[Open File(s)]*** Open interferogram image(s) file(s) with the presence of a gas jet. Image file extensions should preferably be .png or .snp.
However, all image extensions (*.gif*, *.jpg*, *.bmp*, etc) could be used. The path to the opened file is shown in the text box immediately above. If
more than one file has been opened, each file will be analyzed individually, and the average of all results will be presented to the user.
  > **Warning**   
  >  Interferometry Analysis - Gas Jet software only works with grayscale image files. 
  
- ***[Rotate]*** The image rotation in degrees. Positive degrees promote counterclockwise rotation.  

- ***[Original Size]*** Original dimensions of the image file (width,height). 
  > **Note** The interferogram image shown is scaled to screen size (428,342) for users' viewing only. However, all processes to determine the  gas density profile are done with the original dimensions of the image file.

- ***[Interferogram (Ref.)]*** Scaled reference interferogram image.

- ***[Open File]*** Open an undisturbed interferogram image file. Image file extensions should preferably be .png or .snp. However, all image extensions (.gif, .jpg, .bmp, etc) could be used. The path to opened file is shown in text box above. Unlike interferogram LIP files, the algorithm allows the insertion of only one reference file.
  > **Warning**   
  >  Interferometry Analysis - Gas-Jet software only works with grayscale image files. 

- ***[Analyse Data]*** From this command button, the software will apply data processing to generate the accumulated phase-shift map, the radial phase-shift map, and the map of the electron density distribution of the gas.

- ***[Clear]*** Button to clear input and output data.


### Options
- ***[Select Analysis Area]*** Parameters to select the interferogram area to apply the algorithm. The selected area is defined by a rectangle with edges defined by X and Y coordinates (***[Y Coord]*** and ***[X Coord]***). The user that intends to use the whole interferogram needs to uncheck the checkbox  ***[Use select area]***.

- ***[Experimental Parameteres]*** Frame to set the experimental parameters used to obtain the interferogram image. These parameters are:
  - ***[Laser Wavelength]*** ($\lambda$) and ***[Laser bandwidth FHWM]*** in nm;
  - ***[Gas type]*** list box of some types of gases: $H_{2}$, $N_{2}$, $He$ and $Ar$.
  - ***[Polarizability]*** ($\alpha$) in angstrom³. This parameter usually refers to the tendency of matter to acquire an electric dipole moment when subjected to an electric field.
    > **Note**
    > The polarizability value is automatically filled in after selecting the gas type. If the user wants to use gases not yet listed, the gas polarizability value can be entered manually.

- ***[Analysis Parameters]*** Parameters frame to analyze of the interferogram images.
  - ***[Scaling Factor]*** Interferogram image scale in micrometers/pixel.
  - ***[Sigma - Gaussian filter]*** Pixel spread of the gaussian image filter. The initial Sigma depends on the image dimension, but can changed by the
user. 
  - ***[Gaussian Filter Position]*** This parameter is set automatically by the algorithm and this position defines which frequency will be used to apply
the Inverse Fourier Transform and build the phase map of the gas-jet. Both the above parameters are defined in pixels. 
    > **Note** 
    > The algorithm set the frequency that generate a positive phase map. Users can change the filter position.
  
  - ***[Fringes Orientation]*** Definition of the interferogram fringes orientation (vertical or horizontal).
  - ***[Axisymmetric]*** Definition of the axis of symmetry (or axisymmetric) to apply the Inverse Abel Transform. The axisymmetric can be horizontal or vertical.
  - ***[Sigma - Gaussian Blur]*** Spread of the multidimensional gaussian image filter. The standard deviation of the gaussian filter ($\sigma$) defined by the user is equal for all axes..

### Gas-Jet Profile
- ***[Stages]:*** Stages frame allows the visualization of each step of the algorithm.
  - ***[Fourier Transform]*** This image is the 2D Fourier Transform of the interferogram. From this frequency map (Fig. 2.A), the software automatically selects the frequency that generates a positive phase-shift map. The selected frequency is marked by a red line over a pixels line (or column) identifying the ***[Gaussian Filter position]***.  
  > **Note**   
  >  If the ***[Gaussian Filter position]*** is equal to zero, the software will set the new valor automatically.  The user can change this ***[Gaussian Filter position]*** manually.
  - ***[Gaussian Filter]*** This image represents the Gaussian filter map applied to generate the phase map using the selected frequency (Fig. 2.B).

    |<img src = '/Images/Stages1and2.png' width="60%"> |
    |:--:| 
    | *Fig. 2. Example of: (a) 2D frequency domain obtained by the interferogram Fourier Transform with the selected frequency to be filtered; (b) Gaussian filter to be applied on the selected frequency.* |

From the next three steps, users have the option of viewing the 2D maps or 1D curves with standard deviation using the ***[Standard Deviation]*** checkbox.
 
  - ***[Acc. Phase-shift]*** Accumulated phase-shift ($\Delta\phi$) of the plasma (in rad) recovered from the interferograms.
  
    |<img src = '/Images/Stage3.png'>|
    |:--:| 
    | *Fig. 3. Example of: (a) 2D accumulated phase-shift map and (b) 2D standard deviation map; (c) 1D accumulated phase curves and (d) standard deviation of one curve. All phase values are given in* $rad$.|   
    
  - ***[Radial Phase-shift]*** Radial phase-shift ($\Delta\phi_r$) map in $rad/\mu m$ obtained after applying an Inverse Abel Transform to the Accumulated Phase-shift map ($\Delta\phi$).
  
    | ![Phase map](/Images/Stage4.png) |
    |:--:| 
    | *Fig. 4. Example of: (a) 2D radial phase-shift map and (b) 2D standard deviation map; (c) and (d) accuracy between 1D radial phase-shift and normalized phase-shift curves. All phase values are given in* $rad / \mu m$. |  
    
  - ***[Density Profile]*** Molecular density profile ($N$) of the gas jet in $cm^{−3}$ built from the radial phase-shift map ($\Delta\phi_r$) and ***[Experimental Parameters]*** defined by the user: ***[Polarizability]*** ($\alpha$) and ***[Laser Wavelength]*** ($\lambda$).
    
    | ![Gasjet_density](/Images/Stage5.png)|
    |:--:| 
    | *Fig. 5. Example of: (a) 2D gas density map and (b) 2D standard deviation map; (c) 1D gas density curves and (d) standard deviation of one density curve. All density values are given in* $cm^{-3}$. |

- ***[1D Profile]*** This button enables 1D frame (*Fig. 6*) with options for the user to visualize the curves of each select stage for different positions on the chosen symmetrical axis.
- ***[2D Profile]*** This button enables the visualization of each ***[Stage]*** in 2D images.
    
|<img src = '/Images/MainScreen2.png'>|
|:--:| 
| *Fig. 6. Software Main Screen with 1D form enabled and Standard Deviation checkbox selected.* |

- ***[Save Plot]*** This button allows the user to save the visualized plot as an image file (*.png*, *.jpg*, *.bmp*, etc).
- ***[Save Data]*** This button allows the user to save the 2D array that generated the visualized plot as a *.dat* or *.txt* file.
- ***[Colormap dist.]*** With this list box the user can choose between three different color scalings: linear (*Fig 7.a*), quadratic (*Fig 7.a*), or cubic(*Fig 7.c*).

|<img src = '/Images/Colormaps.png' width="100%"> |
|:--:| 
| *Fig. 7. Examples with the colormaps distributions: (a) Linear distribution, (b) Quadratic distribution, (c) Cubic Distribution.* |

## How it works
The detailed description of the algorithm will be presented in a future article. However, the summarized data processing by the software algorithm is described by the flowchart shown in *Fig. 8*:

|<img src = '/Images/Flowchart_gasjet.png'>|
|:--:| 
| *Fig. 8. Scheme of the algorithm data processing.* |

In the scheme of the algorithm data processing (*Fig. 8*): $I$ and $I_0$ are the intensity functions of the bi-dimensional fringes fields obtained from the gas and the reference interferograms, respectively. The hats denotes the Fourier transform of the intensities, $N_gas$ is the gas jet density calculated from its refractive index, $n$, and polarizability, $\alpha$, using the Lorentz-Lorenz relation [[14,15]](#reference).

## Example
In the Example folder of this repository, the user will find the two interferogram images shown in *Fig. 9*. These images were obtained using a Mach-Zehnder-like interferometer, as discussed in [[16]](#reference).

|<img src = '/Example/Interferogram_reference.png' width='40%'> <img src = '/Example/Interferogram_gas.png' width='40%'> |
|:--:| 
| *Fig. 9. Examples of Interferogram images: Reference image (on the left), and gas-jet image (on the rigth).* |

The gas jet in this example was generated by a homemade supersonic nozzle [[17]](#reference) fixed inside a vacuum chamber. During data collection, the nozzle backpressure was maintained at $~50$ $bar$ with nitrogen gas. The interferograms were generated by the second harmonic ($395$ $nm$) of ultrashort pulses from a Ti:sapphire multipass CPA system [[1]](#reference).
All input parameters used for gas-jet characterization are shown in *Fig. 10*.

|<img src = '/Example/MainScreen_example.png'> |
|:--:| 
| *Fig. 10. Software main screen: input parameters and result of the example.* |

## Authors
Interferometry Analysis - Gas-Jet software was developed by researchers of the High-Power Ultrashort Pulse Lasers Group from the Center for Lasers and
Applications (CLA) from the Instituto de Pesquisas Energéticas e Nucleares ([IPEN](https://www.ipen.br/portal_por/portal/default.php)), and of the Instituto Tecnológico de Aeronáutica ([ITA](http://www.ita.br/)).

* Jhonatha Ricardo dos Santos [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7877-0580)
* Armando Valter Felicio Zuffi [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-5705-1499)
* Nilson Dias Vieira Junior [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0003-0092-9357)
* Edison Puig Maldonado [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0002-9462-8151)
* Ricardo Elgul Samad [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7762-8961)

## Acknowledgment
Interferogram Analysis – Gas-Jet was developed to help with the analysis of the gas-jets generated by micro nozzles. This software is part of the work performed by the High-Power Ultrashort Pulse Lasers Group of the CLA/IPEN in the [Extreme Light Laboratory](https://www.unl.edu/diocles/home) (ELL)  at the University of Nebraska - Lincoln (UNL). This partnership was able due to funding provided by the São Paulo Research Foundation (FAPESP) and by the US Department of Energy, through the [LaserNetUS](https://lasernetus.org/) program.

The authors thank all members of the ELL team, especially Prof. Dr. Donald Umstadter director of ELL.

The author Jhonatha Ricardo dos Santos also acknowledges the FAPESP for doctoral fellowship 2021/13737-8. 

## License
Interferogram Analysis - Gas-Jet Profile is licensed under the [MIT license](/LICENSE).

Copyright (c) 2023 Jhonatha Ricardo dos Santos

## Citation
You can find the DOI for the latest version at [Zenodo](https://zenodo.org/badge/latestdoi/601841375).

## Reference
- [1] A. V. F. Zuffi, J. R. dos Santos, E. P. Maldonado, N D. Vieira, and R. E. Samad, "Femtosecond laser-plasma dynamics study by a time-resolved Mach–Zehnder-like interferometer," Appl. Opt. 62, C128-C134 (2023) [DOI: 10.1364/AO.477395](https://doi.org/10.1364/AO.477395).
- [2]  P. Sprangle, B. Hafizi, and J. R. Peñano, “Laser pulse modulation instabilities in plasma channels,” Phys. Rev. E 61, 4381–4393 (2000).[DOI: 10.1103/PhysRevE.61.4381](https://doi.org/10.1103/PhysRevE.61.4381).
- [3] G. Costa, M. P. Anania, F. Bisesto, E. Chiadroni, A. Cianchi, A. Curcio,M. Ferrario, F. Filippi, A. Marocchino, F. Mira, R. Pompili, and A. Zigler,“Characterization of self-injected electron beams from LWFA experiments at SPARC_LAB,” Nucl. Instrum. Methods A 909, 118–122 (2018).[DOI 10.1016/j.nima.2018.02.008](https://doi.org/10.1016/j.nima.2018.02.008).
- [4] G. S. Settles, Schlieren and shadowgraph techniques: visualizing phenomena in transparent media, in Experimental Fluid Mechanics (Springer, 2001), pp. xviii.
- [5] S. Shiraishi, C. Benedetti, A. J. Gonsalves, K. Nakamura, B. H. Shaw, T. Sokollik, J. van Tilborg, C. G. R. Geddes, C. B. Schroeder, C. Toth, E. Esarey, and W. P. Leemans, “Laser red shifting based characterization of wakefield excitation in a laser-plasma accelerator,” Phys. Plasmas 20, 063103 (2013).[DOI 10.1063/1.4810802](https://doi.org/10.1063/1.4810802).
- [6] A. J. Goers, G. A. Hine, L. Feder, B. Miao, F. Salehi, J. K. Wahlstrand, and H. M. Milchberg, “Multi-MeV electron acceleration by Subterawatt laser pulses,” Phys. Rev. Lett. 115, 194802 (2015).[DOI 10.1103/PhysRevLett.115.194802](https://doi.org/10.1103/PhysRevLett.115.194802).
- [7] F. Brandi and L. A. Gizzi, “Optical diagnostics for density measurement in high-quality laser-plasma electron accelerators,” High Power Laser Sci. Eng. 7, e26 (2019).[DOI 10.1017/hpl.2019.11](https://doi.org/10.1017/hpl.2019.11).
- [8] A. K. Arunachalam, “Investigation of laser-plasma interactions at near-critical densities,” Dissertation (University of Jena, 2017)..
- [9] Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). [DOI: 10.1038/s41586-020-2649-2](https://www.nature.com/articles/s41586-020-2649-2). 
- [10] Pauli Virtanen, et. al. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272. [DOI: 10.1038/s41592-019-0686-2](https://www.nature.com/articles/s41592-019-0686-2).
- [11] Gibson, Stephen; Hickstein, Daniel D.; Yurchak, Roman; Ryazanov, Mikhail; Das, Dhrubajyoti; Shih, Gilbert.(2022) PyAbel, PyAbel: v0.9.0, Zenodo,  [DOI: 10.5281/zenodo.7438595](https://doi.org/10.5281/zenodo.7438595).
- [12] Clark, A. (2015). Pillow (PIL Fork) Documentation. readthedocs. Retrieved from [https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf](https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf).
- [13] Stéfan van der Walt, Johannes L. Schönberger, Juan Nunez-Iglesias, François Boulogne, Joshua D. Warner, Neil Yager, Emmanuelle Gouillart,
Tony Yu and the scikit-image contributors. scikit-image: Image processing in Python. PeerJ 2:e453 (2014). [DOI: 10.7717/peerj.453](https://doi.org/10.7717/peerj.453)
- [13] J. D. Hunter, Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering, 9 (3), 90-95 (2007). [
DOI: 10.1109/MCSE.2007.55] (https://ieeexplore.ieee.org/document/4160265)
- [14] H. A. Lorentz, "Über die Beziehungzwischen der Fortpflanzungsgeschwindigkeit des Lichtes derKörperdichte", Ann. Phys. 9, 41-665, (1880). [DOI: 10.1002/andp.18802450406]( https://doi.org/10.1002/andp.18802450406)
- [15] L. Lorenz, "Über die Refractionsconstante", Ann. Phys. 11, 70-103  (1880). [DOI: 10.1002/andp.18802470905](https://doi.org/10.1002/andp.18802470905)
- [16] A. V. F. Zuffi, E. P. Maldonado, N. D. Vieira, and R. E. Samad, “Development of a modified Mach-Zehnder interferometer for time and space density measurements for laser wakefield acceleration,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).
- [17] B. B. Chiomento, A. V. F. Zuffi, N. D. V. Junior, F. B. D. Tabacow, E. P. Maldonado, and R. E. Samad, “Development of dielectric de Laval
nozzles for laser electron acceleration by ultrashort pulses micromachining,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).


