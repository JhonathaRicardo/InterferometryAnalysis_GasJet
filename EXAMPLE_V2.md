# <h1 align = "center"> Using Interferometry Analysis Gas-Jet V.2.0<br> Detailed Example </h1>

## 1. Execute the software
The software can be used in two ways:
A) Use an executable file the executable file available for download [here](https://drive.google.com/file/d/1K_-3wm8TOzxROLyAc22AkRXlRvD3GjjI/view?usp=sharing)
B) Use a Python IDEs and execute de script *IntAnalysis_GasJet_v2.0.py* enabled in this repository.
> NOTE: Make sure all required libraries are installed.   

## 2. Open Files
Use the ***[Open File(s)]*** button to open the interferogram file *Interferogram_gas.png* with the gaseous target and the ***[Open Ref.]*** to open the interferogram *Interferogram_reference.png* used as a reference. Both files are in the *EXAMPLES* folder in this repository. 

> NOTE: The files must have the same dimensions. The user can open N-files with a target and only one reference, or it can open N-files to target and N-files to reference. In this second case, the number of files must be equal. Otherwise, the program will only consider one reference file.

> WARNING: A large number of files will considerably alter processing time. Suggestion: Use less than 10 files!

## 3. Select Area
Use the mouse to select the gaseous target area in the interferogram. After that, use the ***[Y Coord]*** and ***[X Coord]*** for greater precision selection. For this interferogram we used *X Coord = (113, 291)* and *X Coord = (110, 282)* with a *BG Phase Fit* as a median of the 5% smaller phase values.

## 4. Input Parameters

> NOTE: No rotation was necessary. The orientation of the fringes is horizontal.
