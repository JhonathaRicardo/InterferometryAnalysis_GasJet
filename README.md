# <h1 align = "center">Interferometry Analysis - Gas-Jet Profile </h1>
The interferometric technique is an important tool for analysis and diagnosis in astronomy, spectroscopy, metrology, plasma physics, particle physics, and other areas. In Laser Wakefield Acceleration (LWFA) studies, knowing the density distribution of the gas target is crucial to understand the phenomena involved in the particle acceleration process [[1]](#reference).
This Python algorithm was developed to recover the accumulated phase across the neutral gas target as well as estimate its density distribution.
![Logo](https://github.com/JhonathaRicardo/InterferometryAnalysis_GasJetProfile/blob/main/Images/Intro_GasJet.jpg)
![License](https://img.shields.io/badge/license-MIT-green)
![version](https://img.shields.io/badge/version-v.1.0-green)
![status](https://img.shields.io/badge/status-under%20development-green)

## Summary
* [Introduction](#introduction)
* [Installation](#installation)
* [How to use it](#how-to-use-it)
  * [Main Screen](#main-screen)
  * [Interferogram Image](#interferogram-image)
  * [Options](#options)
  * [Gas-Jet Profile](#gas-jet-profile)
* [How it work](#how-it-work)
* [Reference](#reference)
* [Authors](#authors)
* [Acknowledgment](#acknowledgment)
 * [License](#license)
* [Citation](#citation)

## Introduction
In recent decades, the continuous development of compact particle accelerators based on wakefield laser acceleration (LWFA) has promoted contributions to fundamental and applied research [[1, 2]](#reference), including possible future use in proton therapy and hadron therapy [[3, 4]](#reference), and in the production of radioisotopes for nuclear medicine [[5, 6]](#reference). In this approach, ultrafast high-intensity laser pulses focused on a gas target to create a plasma wave with longitudinal electric fields able to accelerate electrons [[7, 8]](#reference). The advances in high-peak-power tabletop lasers in the last years and the high longitudinal electric fields (up to ∼1 TV/m) supported by plasma waves attracted attention to LWFA as a compact alternative for RF conventional accelerators [[9]](#reference).

Many groups around the world have sought advances in the field of LWFA and other plasma acceleration schemes from institutions in North America, Europe, and Asia [[2]](#reference). In Latin America, our research group pioneered the implementation of a laser-plasma accelerator at the Institute of Energy and Nuclear Research (IPEN) [[10]](#reference).Our main objective is to produce electron beams with energy up to tens of megaelectron volts per LWFA. By bremsstrahlung, Electron beams with those energies are able to produce $\gamma$ radiation to induce a $^{100} Mo(\gamma , n) ^{99} Mo$ photonuclear reaction as a future application [[11]](#reference). 

Currently, we are focusing efforts on several developments required for a LWFA installation, such as computational simulation support [[12–14]](#reference), a source of high-peak-power laser pulses [[15]](#reference), proper gaseous and plasma target creation [[16, 17]](#reference), and development and implementation of diagnostic tools to assist and monitor the experiments [[18, 19]](#reference). The development of diagnostic tools is very important for a better understanding of the laser-plasma interaction [[20]](#reference). Diagnostic efficiency is crucial, as instabilities in both targets and laser pulses can result in low reproducibility of LWFA processes and impair the quality of accelerated electron beams [[21]](#reference).
Among the various non-perturbing optical methods that can be used to diagnose the gaseous target [[22-25]](#reference), interferometry is a very accurate technique capable of quantifying very small optical path differences and therefore suitable for measuring density variations in LWFA targets [[26, 27]](#reference). 

The [Interferometry Analysis - Gas-Jet Profile] software and [Interferometry Analysis - LIP Profile] software were developed due to this need for a new diagnostic tool to aid in the characterization of the supersonic jet of gas, quickly and reliably. Both softwares were developed by our research group as part of the work to implementation of a laser-plasma accelerator at the Nuclear and Energy Research Institute (IPEN).

## Installation
Interferometry Analysis - Gas-Jet software was developed in Python 3.11 and the use of this algorithm requires the installation of the following packages: [NumPy](https://numpy.org/) [[2]](#reference), [Scipy](https://scipy.org/) [[3]](#reference) and [PyAbel](https://pyabel.readthedocs.io/en/latest/index.html) [[4]](#reference) for data processing, [Pillow](https://pypi.org/project/Pillow/) [[5]](#reference) to manipulate interferogram images, [Matplotlib](https://matplotlib.org/stable/index.html) [[6]](#reference) to plot results, and
[PySimpleGui](https://www.pysimplegui.org/en/latest/) to create the users template.

The second way to use this software is through the executable file. The users can create a single .exe file using the [pyinstaller](https://pyinstaller.org/en/stable/) package trought the follow terminal command:

<code>   pyinstaller --onefile -w IntAnalysis_GasJetProfile.py                </code>

## How to use it
Interferometry Analysis - The Gas-Jet software has a graphical interface developed with PysimpleGUI. This interface assists users and facilitates their applications.
In this section, we provide users with a simple review of the software's functions and how to use them.
### Main Screen

| ![Template](/Images/Template_image.jpg) |
|:--:| 
| *Fig.1 - Software Template* |

### Interferogram Image
- ***[Interferogram (Gas-Jet)]*** Scaled gas-jet interferogram image.

- ***[Open File(s)]*** Open interferogram image(s) file(s) with the presence of gas jet. Image file extensions should preferably be .png or .snp. However, all image extensions (.gif, .jpg, .bmp, etc) could be used. The path to opened file is shown in text box above. If more than one file has been opened, each file will be analyzed individually and the average of all results will be presented to the user.
  > **Warning**   
  >  Interferometry Analysis - Gas-Jet software only works with grayscale image files. 

- ***[Interferogram (Ref.)]*** Scaled reference interferogram image.

- ***[Open File]*** Open an undisturbed interferogram image file. Image file extensions should preferably be .png or .snp. However, all image extensions (.gif, .jpg, .bmp, etc) could be used. The path to opened file is shown in text box above. Unlike interferogram gas jet files, the algorithm allows the insertion of only one reference file.
  > **Warning**   
  >  Interferometry Analysis - Gas-Jet software only works with grayscale image files. 
 
- ***[Rotate]*** The image rotates in degrees. Positive degrees promote counterclockwise rotation.  

- ***[Image Scale]*** The interferogram image shown is scaled to screen size (428,342) for users' viewing only. However, all processes to determine the gas jet density profile are done with the original dimensions of the image file.

- ***[Analyse Data]*** From this command button, the software will apply data processing to generate accumulated phase, inverse Abel transforms, and gas jet density profile.

- ***[Exit]*** Close software.

### Options
- ***[Select Analysis Area]*** From the parameters in this form, the user can select the interferogram area to apply the algorithm to determine the gas-jet density profile. The selected area is defined by a rectangle with edges defined by X and Y coordinates (***[Y Coord]*** and ***[X Coord]***).
The user that intends to use the whole interferogram figure needs to uncheck the checkbox ***[Use select area]***.

- ***[Experimental Parameteres]*** In this form, the user set the experimental parameters used to obtain the interferogram image. These parameters are:
  - ***[Laser Wavelength]*** and ***[uncertainty Laser Wavelength]*** in nm;
  - Gas ***[Polarizability]*** in angstrom³. This parameter usually refers to the tendency of matter to acquire an electric dipole moment when subjected to an electric field.

- ***[Analysis Parameters]*** This form contains the parameters for analysis of the interferogram images:
  - ***[Scaling Factor]*** of an interferogram image in pixels/micrometers;
  - ***[Sigma - Gaussian Blur]*** is a multidimensional gaussian image filter. The standard deviation of the gaussian filter (Sigma) defined by the user is equal for all axes.
  - ***[ Gaussian Filter Position]*** this parameter is different from Gaussian Blur. This parameter is set automatically by the algorithm and this position defines which frequency will be used to apply the Inverse Fourier Transform and build the phase map of the gas-jet.
Both the above parameters are defined in pixels. 
    > **Note** 
    > The algorithm set the frequency that defines a positive phase map.  But, users can change the filter position.
  
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

From the next three steps, users have the option of viewing the average values of the maps in 2D (using ***[2D Profile]*** button) or the 1D profile of the maps in different positions on a symmetrical axis (using ***[1D Profile]*** button). 2D standard deviation maps or 1D curves can be viewed using the ***[Standard Deviation]*** checkbox.
  
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
    
    
## How it work
The interferogram analysis software algorithm works according to the flowchart below:

### Accumulated Phase
The Accumulated Phase map (or accumulated phase shift map) is obtained from the shifts of the speckle fields from two interferogram images. The first is the interferogram image with fringes disturbed due to the presence of gas and the second is a backgroung image (or reference) with undisturbed fringes. According to the flowchart, apply 2D Fourier transforms on both interferograms by transporting them in the frequency domain. Applying a Gaussian filter over the region containing the phase shift information [[7]](#reference) and inverting the Fourier transform over two frequency domain maps. Finally, we obtain the accumulated (or integrated) phase-shift map $&Delta;&phi;_{z}$ [[7, 8]](#reference) along the beam propagation direction (z direction) by the following equation:

$$ 
\begin{equation}
\Delta\phi_{z} =  tan^{-1}\left({\phi_{gas} - \phi_{ref}}\right)
\tag{1}
\end{equation}
$$

where  $\phi_{gas}$  and  $\phi_{ref}$   is, respectively, the gas-jet and background phase map.

#### Standard Deviation of Accumulated Phase

According to M. Lehmann [[9]](#reference), for two well-resolved speckle fields (background and perturbed by gas) the phase error is determined by the probability distributions of the intensities and phase derivatives of two speckle fields. Considering that each field has a Gaussian distribution of speckle intensities and since the measured phase is the difference between two speckle phases, its error also follows a Gaussian probability distribution, with standard deviation given by:

$$ 
\begin{equation}
\sigma_{\Delta\phi}(\Delta x, I_{1} ,I_{2}) = {\Delta x \over 2}{\pi \over \beta} \left\lbrack{ \overline{I} (I_{1} + I_{2}) \over 2 I_{1} I_{2}}\right\rbrack^{1/2}
\tag{2}
\end{equation}
$$

where $I_{1}$ and $I_{2}$ are the intensity distribution of both speckle fields (with mean intensity $\overline{I}$ each), $&Delta;x$ is the displacement of the image in the direction perpendicular to the direction axis of the fringes, and $\beta$ is the speckle size. 

### Inverse Abel Transform

As mentioned above, $\Delta\phi_{z}$ the integrated phase map along the laser beam propagation direction (z direction).
Assuming an axisymmetric gas-jet, the integrated information along $z$ is sufficient to reconstruct the radial information using inversion $\Delta\phi_{r}$ such as the Abel inversion method [[8]](#reference).

$$ \Delta\phi_{r} = - {1 \over \pi} \int_{r}^{\infty} {d (\Delta\phi_{z}) \over dz} {dz \over \sqrt {z² - r²}} \tag{3}$$ 

In this software, the phase map is determined from the application of the PyAbel [[4]](#reference) algorithm on the accumulated phase map.
PyAbel is a Python package that provides functions for the forward and inverse Abel transforms. The inverse Abel transform takes a 2D projection and reconstructs a slice of the cylindrically symmetric 3D distribution, which makes this function an important tool in analyzing the projections of angle-resolved, plasma plumes, flames, solar occultation [[4]](#reference), and gas-jets.

PyAbel provides efficient implementations of several Abel transform algorithms [[10]](#reference). In this software, the chosen method was the deconvolution algorithm Dash Onion Peeling because it is simple and computationally very efficient [[4]](#reference). According to Dash [[11]](#reference), this method requires less smoothing than other methods.

#### Standard Deviation of Inverse Abel Transform
The accuracy of applying the inverse Abel transform is associated with the standard deviation generated by a convolution between the phase map and the integrated phase map.

### Density Profile
The gas molecular density $\rho_{gas}$ is obtained through the simple expression of the Lorentz-Lorenz relation [[12, 13]](#reference):

$$ \rho_{gas} = {3 \over 4\pi\alpha} {(n^2-1) \over (n^3+2)} \tag{4}$$

where $\alpha$ is the molecular polarizability of gas and $n$ is the spatial refractive index distribution. This 2D refractive index map is determined from the phase-shift map $\Delta\phi_{r}$:

$$ n = 1 + {\Delta\phi_{r} \lambda \over 2\pi} \tag{4}$$



#### Standard Deviation of Density
The accuracy of the gas density measurement depends on the accuracies of the phase-shift measurement ($\sigma_{\Delta\phi_{z}}$), and the numerical accuracy of the Abel inversion [[14]](#reference) ($\sigma_{abel}$). This way, the standard deviation of phase-shift map $\sigma_{\Delta\phi_{r}}$ can be write as ($5$):

$$  \sigma_{\Delta\phi_{r}} = \sqrt{\left({\sigma_{\Delta\phi_{z}}}^2 + {\sigma_{abel}}^2\right)} \tag{5}$$

So, the standard deviation of gas density $\sigma_{\rho}$ is given by ($6$):

$$ 
\begin{equation}
\sigma_{\rho} = \sqrt{\left({\partial\rho \over \partial\Delta\phi_{r}}\right)^2 \left({\sigma_{\Delta\phi_{r}}}\right)^2 + 
\left({\partial\rho \over \partial\lambda}\right)^2 \left({\sigma_{\lambda}}\right)^2}
\tag{6}
\end{equation}
$$

where $\sigma_{\lambda}$ is user-defined experimental parameter.
  > **Note**
  > The contribution of this parameter to the density standard deviation is usually very small. 
  
## Reference
- [1] F. Albert, “Laser wakefield accelerators: next-generation light sources,” Opt. Photonics News 29(1), 42–49 (2018). [DOI: 10.1364/OPN.29.1.000042](https://doi.org/10.1364/OPN.29.1.000042).
- [2] F. Albert, M. E. Couprie, A. Debus, et al., “2020 roadmap on plasma accelerators,” New J. Phys. 23, 031101 (2021). [DOI: 10.1088/1367-2630/abcc62](https://iopscience.iop.org/article/10.1088/1367-2630/abcc62).
- [3] K. Ledingham, P. Bolton, N. Shikazono, and C. M. Ma, “Towards laser driven Hadron cancer radiotherapy: a review of progress,” Appl. Sci. 4, 402–443 (2014). [DOI: 10.3390/app4030402](https://doi.org/10.3390/app4030402).
- [4] A. Giulietti, Laser-Driven Particle Acceleration Towards Radiobiology and Medicine, Biological and Medical Physics, Biomedical Engineering (Springer International Publishing, 2016).
- [5] K. Nemoto, A. Maksimchuk, S. Banerjee, K. Flippo, G. Mourou, D. Umstadter, and V. Y. Bychenkov, “Laser-triggered ion acceleration and table top isotope production,” Appl. Phys. Lett. 78, 595–597 (2001). [DOI: 10.1063/1.1343845](https://doi.org/10.1063/1.1343845).
- [6] I. Spencer, K. W. D. Ledingham, R. P. Singhal, T. McCanny, P. McKenna, E. L. Clark, K. Krushelnick, M. Zepf, F. N. Beg, M. Tatarakis, A. E. Dangor, P. A. Norreys, R. J. Clarke, R. M. Allott, and I. N. Ross, “Laser generation of proton beams for the production of short-lived positron emitting radioisotopes,” Nucl. Instrum. Methods B 183, 449–458 (2001). [DOI: 10.1016/S0168-583X(01)00771-6](https://doi.org/10.1016/S0168-583X(01)00771-6).
- [7] E. Esarey, C. B. Schroeder, and W. P. Leemans, “Physics of laser-driven plasma-based electron accelerators,” Rev. Mod. Phys. 81, 1229–1285 (2009). [DOI: 10.1103/RevModPhys.81.1229](https://doi.org/10.1103/RevModPhys.81.1229).
- [8] T. Tajima and J. M. Dawson, “Laser electron-accelerator,” Phys. Rev. Lett. 43, 267–270 (1979).[DOI: 10.1103/PhysRevLett.43.267](https://doi.org/10.1103/PhysRevLett.43.267).
- [9] S. M. Hooker, “Developments in laser-driven plasma accelerators,” Nat. Photonics 7, 775–782 (2013).[DOI: 10.1038/nphoton.2013.234](https://doi.org/10.1038/nphoton.2013.234).
- [10] N. D. Vieira, R. E. Samad, and E. P. Maldonado, “Compact laser accelerators towards medical applications—perspectives for a Brazilian Program,” in SBFoton International Optics and Photonics Conference (IEEE, 2019).
- [11] N. D. Vieira, E. P. Maldonado, A. Bonatto, R. P. Nunes, S. Banerjee, F. A. Genezini, M. Moralles, A. V. F. Zuffi, and R. E. Samad, “Laser wake-
field electron accelerator: possible use for radioisotope production,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).
- [12] E. P. Maldonado, R. E. Samad, A. Bonatto, R. P. Nunes, S. Banerjee, and N. D. Vieira, “Study of quasimonoenergetic electron bunch generation in self-modulated laser wakefield acceleration using TW or sub-TW ultrashort laser pulses,” AIP Adv. 11, 065116 (2021).[DOI: 10.1063/5.0052831](https://doi.org/10.1063/5.0052831).
- [13] E. P. Maldonado, R. E. Samad, A. Bonatto, R. P. Nunes, S. Banerjee, and N. D. Vieira, “Electron beam properties in self-modulated laser wakefield acceleration using TW and sub-TW pulses,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).
- [14] E. P. Maldonado, R. E. Samad, A. V. F. Zuffi, F. B. D. Tabacow, and N. D. Vieira, “Self-modulated laser-plasma acceleration in a H2 gas target, simulated in a spectral particle-in-cell algorithm: wakefield and electron bunch properties,” in SBFoton International Optics and Photonics Conference (IEEE, 2019).
- [15] R. E. Samad, E. P. Maldonado, W. De Rossi, and N. D. V. Junior, “High intensity ultrashort laser pulses and their applications at IPEN,”
in SBFoton International Optics and Photonics Conference (IEEE,2021).
- [16] B. B. Chiomento, A. V. F. Zuffi, N. D. V. Junior, F. B. D. Tabacow, E. P. Maldonado, and R. E. Samad, “Development of dielectric de Laval
nozzles for laser electron acceleration by ultrashort pulses micromachining,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).
- [17] F. B. D. Tabacow, A. V. F. Zuffi, E. P. Maldonado, R. E. Samad, and N. D. Vieira, “Theoretical and experimental study of supersonic gas jet targets for laser wakefield acceleration,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).
- [18] A. V. F. Zuffi, E. P. Maldonado, N. D. Vieira, and R. E. Samad, “Development of a modified Mach-Zehnder interferometer for time
and space density measurements for laser wakefield acceleration,” in SBFoton International Optics and Photonics Conference (IEEE, 2021).
- [19] R. E. Samad, A. V. F. Zuffi, E. P. Maldonado, and N. D. Vieira, “Development and optical characterization of supersonic gas targets for high-intensity laser plasma studies,” in SBFoton International Optics and Photonics Conference (IEEE, 2018).
- [20] A. V. F. Zuffi, J. R. dos Santos, E. P. Maldonado, N D. Vieira, and R. E. Samad, "Femtosecond laser-plasma dynamics study by a time-resolved Mach–Zehnder-like interferometer," Appl. Opt. 62, C128-C134 (2023) [DOI: 10.1364/AO.477395](https://doi.org/10.1364/AO.477395).
- [21]  P. Sprangle, B. Hafizi, and J. R. Peñano, “Laser pulse modulation instabilities in plasma channels,” Phys. Rev. E 61, 4381–4393 (2000).[DOI: 10.1103/PhysRevE.61.4381](https://doi.org/10.1103/PhysRevE.61.4381).
- [22] G. Costa, M. P. Anania, F. Bisesto, E. Chiadroni, A. Cianchi, A. Curcio,M. Ferrario, F. Filippi, A. Marocchino, F. Mira, R. Pompili, and A. Zigler,“Characterization of self-injected electron beams from LWFA experiments at SPARC_LAB,” Nucl. Instrum. Methods A 909, 118–122 (2018).[DOI 10.1016/j.nima.2018.02.008](https://doi.org/10.1016/j.nima.2018.02.008).
- [23] G. S. Settles, Schlieren and shadowgraph techniques: visualizing phenomena in transparent media, in Experimental Fluid Mechanics (Springer, 2001), pp. xviii.
- [24] S. Shiraishi, C. Benedetti, A. J. Gonsalves, K. Nakamura, B. H. Shaw, T. Sokollik, J. van Tilborg, C. G. R. Geddes, C. B. Schroeder, C. Toth, E. Esarey, and W. P. Leemans, “Laser red shifting based characterization of wakefield excitation in a laser-plasma accelerator,” Phys. Plasmas 20, 063103 (2013).[DOI 10.1063/1.4810802](https://doi.org/10.1063/1.4810802).
- [25] A. J. Goers, G. A. Hine, L. Feder, B. Miao, F. Salehi, J. K. Wahlstrand, and H. M. Milchberg, “Multi-MeV electron acceleration by Subterawatt laser pulses,” Phys. Rev. Lett. 115, 194802 (2015).[DOI 10.1103/PhysRevLett.115.194802](https://doi.org/10.1103/PhysRevLett.115.194802).
- [26] F. Brandi and L. A. Gizzi, “Optical diagnostics for density measurement in high-quality laser-plasma electron accelerators,” High Power Laser Sci. Eng. 7, e26 (2019).[DOI 10.1017/hpl.2019.11](https://doi.org/10.1017/hpl.2019.11).
- [27] A. K. Arunachalam, “Investigation of laser-plasma interactions at near-critical densities,” Dissertation (University of Jena, 2017)..


- [1] Hariharan, P. (2007) Basics of Interferometry. 2nd Edition, Elsevier, Amsterdam.[DOI: 10.1016/B978-0-12-373589-8.X5000-7](https://doi.org/10.1016/B978-0-12-373589-8.X5000-7).
- [2] Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). [DOI: 10.1038/s41586-020-2649-2](https://www.nature.com/articles/s41586-020-2649-2). 
- [3] Pauli Virtanen, et. al. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272. [DOI: 10.1038/s41592-019-0686-2](https://www.nature.com/articles/s41592-019-0686-2).
- [4] Gibson, Stephen; Hickstein, Daniel D.; Yurchak, Roman; Ryazanov, Mikhail; Das, Dhrubajyoti; Shih, Gilbert.(2022) PyAbel, PyAbel: v0.9.0, Zenodo,  [DOI: 10.5281/zenodo.7438595](https://doi.org/10.5281/zenodo.7438595).
- [5] Clark, A. (2015). Pillow (PIL Fork) Documentation. readthedocs. Retrieved from [https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf](https://buildmedia.readthedocs.org/media/pdf/pillow/latest/pillow.pdf).
- [6] J. D. Hunter, Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering, 9 (3), 90-95 (2007). [
DOI: 10.1109/MCSE.2007.55] (https://ieeexplore.ieee.org/document/4160265)
- [7] J. P. Couperus, A. Kohler, T. A. W. Wolterink, A. Jochmann, O. Zarini, H. M. J. Bastiaens, K. J. Boller, A. Irman, and U. Schramm, Nucl Instrum Meth A 830, 504-509 (2016).[DOI: 10.1016/j.nima.2016.02.099](https://doi.org/10.1016/j.nima.2016.02.099).
- [8] V. Malka, C. Coulaud, J. P. Geindre, V. Lopez, Z. Najmudin, D. Neely, and F. Amiranoff, Rev. Sci. Instrum. 71, 2329-2333 (2000). [DOI: 10.1063/1.1150619](https://doi.org/10.1063/1.1150619)
- [9] Mathias Lehmann, "Decorrelation-induced phase errors in phase-shifting speckle interferometry," Appl. Opt. 36, 3657-3667 (1997). [DOI: 10.1364/AO.36.003657](https://doi.org/10.1364/AO.36.003657).
- [10] Daniel D. Hickstein, Stephen T. Gibson, Roman Yurchak, Dhrubajyoti D. Das, Mikhail Ryazanov. A direct comparison of high-speed methods for the numerical Abel transform. Rev. Sci. Instrum., 90, 065115, 2019. [DOI: 10.1063/1.5092635](https://doi.org/10.1063/1.5092635).
- [11] C. J. Dasch, “One-dimensional tomography: a comparison of Abel, onion-peeling, and filtered backprojection methods”, Appl. Opt. 31, 1146–1152 (1992). [DOI: 10.1364/AO.31.001146](https://doi.org/10.1364/AO.31.001146).
- [12] H. A. Lorentz, "Über die Beziehungzwischen der Fortpflanzungsgeschwindigkeit des Lichtes derKörperdichte", Ann. Phys. 9, 41-665, (1880). [DOI: 10.1002/andp.18802450406]( https://doi.org/10.1002/andp.18802450406)
- [13] L. Lorenz, "Über die Refractionsconstante", Ann. Phys. 11, 70-103  (1880). [DOI: 10.1002/andp.18802470905](https://doi.org/10.1002/andp.18802470905)
- [14] A. Saville, M. (2022). 2D Relative Phase Reconstruction in Plasma Diagnostics. Optical Interferometry - A Multidisciplinary Technique in Science and Engineering. [DOI: 10.5772/intechopen.104748](https://www.intechopen.com/chapters/81777).


## Authors
Interferometry Analysis - Gas-Jet software was developed by researchs of the High Power Ultrashort Pulse Lasers Group of the Center for Lasers and Applications (CLA) at Instituto de Pesquisas Energéticas e Nucleares ([IPEN](https://www.ipen.br/portal_por/portal/default.php)).
* Jhonatha Ricardo dos Santos [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7877-0580)
* Armando Valter Felicio Zuffi [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-5705-1499)
* Ricardo Elgul Samad [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0001-7762-8961)
* Nilson Dias Vieira Junior [![logo_ORCID](/Images/logo_ORCID.png)](https://orcid.org/0000-0003-0092-9357)

## Acknowledgment
Interferogram Analysis - Gas-Jet software was developed to help with the analyze the density of supersonic gas-jet generated by micro-nozzles. This software is part of the work realized by the High Power Ultrashort Pulse Lasers Group of the CLA/IPEN in the [Extreme Light Laboratory](https://www.unl.edu/diocles/home) (ELL)  at the University of Nebraska - Lincoln (UNL). This partnership was able due to funding provided by the São Paulo Research Foundation (FAPESP) and by the US Department of Energy, through the [LaserNetUS](https://lasernetus.org/) program.

The author Jhonatha Ricardo dos Santos also acknowledges the FAPESP for doctoral fellowship 2017/13737-8. 

## License
Interferogram Analysis - Gas-Jet Profile is licensed under the [MIT license](/LICENSE).

Copyright (c) 2023 Jhonatha Ricardo dos Santos

## Citation

