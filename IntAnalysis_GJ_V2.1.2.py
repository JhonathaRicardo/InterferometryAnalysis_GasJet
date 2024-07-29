# Software: Interferometry Analysis - Gas Jet - (Version 2.1.2)
# Authors: Jhonatha Ricardo dos Santos, Armando Zuffi, Ricardo Edgul Samad, Nilson Dias Vieira Junior
# Python 3.11
# Last update: 2024_07_29

# LYBRARIES
# The Python Standard Library
# PyAbel/PyAbel:v0.9.0rc1 from https://doi.org/10.5281/zenodo.7401589.svg
# PySimpleGUI from pysimplegui.org
# Matplotlib from matplotlib.org
# Scipy from scipy.org
# Scikit-image from  https://doi.org/10.7717/peerj.453
# Pillow (PIL Fork) 9.3.0 from pypi.org/project/Pillow
import abel
import PySimpleGUI as sg
import os
import io
import math
import numpy as np
import shutil
import tempfile
import matplotlib
import matplotlib.pyplot as plt
import warnings

from io import BytesIO
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.ndimage import gaussian_filter, rotate
from scipy.signal import peak_widths, find_peaks
from scipy import sparse
from scipy.sparse.linalg import spsolve
from PIL import Image, ImageDraw, UnidentifiedImageError
from skimage.restoration import unwrap_phase
from skimage.registration import phase_cross_correlation

# Matplotlib Tk style
matplotlib.use('TkAgg')
warnings.filterwarnings("ignore")
#Colormaps
custom1_cmap = LinearSegmentedColormap.from_list('default',
                                             [(0,    'white'),
                                              (0.25, 'blue'),
                                              (0.40, 'green'),
                                              (0.55, 'lime'),
                                              (0.70, 'yellow'),
                                              (0.85, 'orange'),
                                              (1, 'red')], N=256)
matplotlib.colormaps.register(cmap=custom1_cmap)

# Font and theme of PysimpleGUI
AppFont = 'Arial 16 bold'
sg.theme('DarkGrey4')

# Image files types
file_types = [("SNP (*.snp)", "*.snp"), ("PNG (*.png)", "*.png"), ("All files (*.*)", "*.*")]
# Temp files
tmp_file = tempfile.NamedTemporaryFile(suffix=".png").name
tmp_file2 = tempfile.NamedTemporaryFile(suffix=".png").name
tmp_file_plot = 'temp_plot_abel.png'

# Image files types
file_types = [("SNP (*.snp)", "*.snp"), ("PNG (*.png)", "*.png"), ("All files (*.*)", "*.*")]
# Temp files
tmp_file = tempfile.NamedTemporaryFile(suffix=".png").name
tmp_file2 = tempfile.NamedTemporaryFile(suffix=".png").name
tmp_file_plot = 'temp_plot_abel.png'

# INITIALPARAMETERS
# Image paths
path1 = ''
path2 = ''
# Physics Parametres
lambda0 = '395'  # nm
unc_lambda0 = '0'
factor = '1.000'  # factor um/pixel
unc_factor = '0.001'
polargas = '1.710'  # for N2 gas in A^3
sigma_gfilter = '0'  # sigma of gauss function
sigma_gblur = '5'  # sigma of gaussian blur
centerfilter = '0'  # Center of the gaussian filter application
axis_pos = '0'  # Axisymmetrical position pixel
base_ref = '1'
# Image parameters
h_prof = -1.0  # heigth null
rotate_degree = 0.0  # angle to image rotation
# Initial values to cut image
begin_x = '100'
begin_y = '100'
end_x = '300'
end_y = '300'
# Initial values of heigths for 1D analysis
pos1 = '10'
pos2 = '20'
pos3 = '30'

# Images Dimensions
width, height = size = 428, 342  # Scale image - interferogram
width2, height2 = size2 = 214, 172  # Scale image - Ref
width3, height3 = size3 = 428, 342  # Scale image - Result
# Min and Max values of Interferogram Image
minvalue_x, maxvalue_x, minvalue_y, maxvalue_y = 0, 428, 0, 342
# Frame 1D visible
visible_f1d = False


#################################################################################
# FUNCTIONS
################################################################################
# GET BINARY DATA
def getBinaryData(filename):
    '''
    path file name to binary value
    :param filename: path of file
    :return: binary value
    '''
    binary_values = []
    with open(filename, 'rb') as f:
        data = f.read(1)
        while data != b'':
            binary_values.append(ord(data))
            data = f.read(1)
        return binary_values


# DRAW FIGURE FROM FILES
def draw_figure(canvas, figure):
    '''
    Drawing rectangle figure on canvas
    :param canvas: image canvas
    :param figure: original interferogram
    :return: rectangle drawn on figure
    '''
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


# GET VALUES FROM INPUTBOX
def get_value(key, values):
    '''
    convert string labels to float values.
    :param key: labels
    :param values: labels value
    :return: float of label values
    '''
    value = values[key]
    return float(value)


# GET IMAGE FILE FROM PATH FILE
def image_to_data(im):
    '''
    convert image to data image
    :param im: image
    :return: data of image
    '''
    with BytesIO() as output:
        im.save(output, format="PNG")
        data = output.getvalue()
    return data


# DRAW RECTANGLE ON INTERFEROMETER FIGURE
def apply_drawing(values, window):
    '''
    :param values: x and y labels of rectangle
    :param window: main window
    :return: rectangle drown on temp image file
    '''
    image_file = values["file1"]
    begin_x = get_value("-BEGIN_X-", values)
    begin_y = get_value("-BEGIN_Y-", values)
    end_x = get_value("-END_X-", values)
    end_y = get_value("-END_Y-", values)
    rotate_degree = get_value("-DEGREE-", values)

    if begin_x > end_x:
        begin_x = get_value("-END_X-", values)
        end_x = get_value("-BEGIN_X-", values)
    if begin_y > end_y:
        begin_y = get_value("-END_Y-", values)
        end_y = get_value("-BEGIN_Y-", values)

    if os.path.exists(image_file):
        shutil.copy(image_file, tmp_file)
        imagetmp = Image.open(tmp_file)
        imagetmp = imagetmp.resize(size)
        imagetmp = imagetmp.rotate(rotate_degree, resample=Image.Resampling.BICUBIC)
        draw = ImageDraw.Draw(imagetmp)
        draw.rectangle((begin_x, begin_y, end_x, end_y), width=2, outline='white')  ##DCDCDC
        imagetmp.save(tmp_file)
        bio = io.BytesIO()
        imagetmp.save(bio, format='PNG')
        window["image1"].update(data=bio.getvalue(), size=size)

def rotate_img(values, window):
    '''
    :param values: x and y labels of rectangle
    :param window: main window
    :return: rectangle drown on temp image file
    '''
    image_file = values["file1"]
    rotate_degree = get_value("-DEGREE-", values)

    if os.path.exists(image_file):
        shutil.copy(image_file, tmp_file)
        imagetmp = Image.open(tmp_file)
        imagetmp = imagetmp.resize(size)
        imagetmp = imagetmp.rotate(rotate_degree, resample=Image.Resampling.BICUBIC)
        imagetmp.save(tmp_file)
        bio = io.BytesIO()
        imagetmp.save(bio, format='PNG')
        window["image1"].update(data=bio.getvalue(), size=size)

# CREATING MEAN MAPS/ARRAY AND STD ARRAY
def mean_maps(data):
    '''
    2D Array mean
    :param n: group of 2D arrays
    :return: 2D array
    '''
    mean_data = data[0] / len(data)
    for i in range(1, len(data)):
        mean_data = mean_data + data[i] / len(data)
    return mean_data


def std_maps(data, mean_data):
    '''
    standard deviation of 2D Array maps
    :param n: group of 2D arrays and mean 2D array
    :return: std of 2D array
    '''
    desv = (data[0] - mean_data) * (data[0] - mean_data)
    for i in range(1, len(data)):
        desv = desv + (data[i] - mean_data) * (data[i] - mean_data)

    return np.sqrt(desv / len(data))


# CREATING FRINGES WIDTHS
def fringes_width(data1):
    '''
    Calculate 2D array shifts and widths fringes distribution
    :param n: 2D array, 2D array of ref. image.
    :return: mean fringe width
    '''
    # array with fringes width
    f_width = np.zeros(np.shape(data1))
    nl, nr = np.shape(data1)
    for l in range(0, nl):
        line1 = (data1[l, :])
        ypeaks1, _ = find_peaks(line1)
        y = np.diff(ypeaks1)
        x = np.linspace(0, nr, len(y))
        x_interp = np.linspace(0, nr, nr)
        f_width[l, :] = np.interp(x_interp, x, y)

    return f_width

def baseline2D(data, base_ref):
    '''
    path file name to binary value
    :param 2D array and base ref (%)
    :return: 2D array BG map and 2D array std of BG map
    '''
    data_sort = np.sort(data.flatten())
    bg_map = np.mean(data_sort[0:int(len(data_sort)* base_ref)]) * np.ones(np.shape(data))
    bg_std = np.std(data_sort[0:int(len(data_sort)* base_ref)]) * np.ones(np.shape(data))

    return bg_map, bg_std

def als(data, lam=5e2, p=0.1, itermax=10):
    r"""
    Implements an Asymmetric Least Squares Smoothing
    baseline correction algorithm (P. Eilers, H. Boelens 2005)
    Inputs:
        y:
            input data (i.e. chromatogram of spectrum)
        lam:
            parameter that can be adjusted by user. The larger lambda is,
            the smoother the resulting background, z
        p:
            wheighting deviations. 0.5 = symmetric, <0.5: negative
            deviations are stronger suppressed
        itermax:
            number of iterations to perform
    Output:
        the fitted background vector

    """
    L = len(data)
#  D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    D = sparse.eye(L, format='csc')
    D = D[1:] - D[:-1]  # numpy.diff( ,2) does not work with sparse matrix. This is a workaround.
    D = D[1:] - D[:-1]
    D = D.T
    w = np.ones(L)
    for i in range(itermax):
        W = sparse.diags(w, 0, shape=(L, L))
        Z = W + lam * D.dot(D.T)
        z = spsolve(Z, w * data)
        w = p * (data > z) + (1 - p) * (data < z)
    return z

def func_gfilter(data, centerfilter, f_range, sigma_gfilter):
    '''
    2D gaussian filter
    :param 2D array, filter position, filter range and sigma of Gfilter
    :return: 2D arrayfilter and sigma of Gfilter
    '''
    gfilter = np.zeros(np.shape(data))
    nlines, nrows = np.shape(data)
    X = np.arange(0, nrows, 1)
    Y = np.arange(0, nlines, 1)
    # Creating Filter for Horizontal/vertical fringes orientation
    if values['-combofringe-'] == 'horizontal':
        if sigma_gfilter == 0:
            # sigma filter is a func of image dimensions and f_rqnge
            sigma_gfilter = (2 * f_range)
            window['-sigma_gfilter-'].update(str(sigma_gfilter))
        for i in X:
            gfilter[:, i] = np.exp(-np.square(Y - centerfilter) / (2 * np.square(sigma_gfilter)))

    elif values['-combofringe-'] == 'vertical':
        if sigma_gfilter == 0:
            # sigma filter is a func of image dimensions and f_rqnge
            sigma_gfilter = 2 * f_range
            window['-sigma_gfilter-'].update(str(sigma_gfilter))
        for i in Y:
            gfilter[i, :] = np.exp(-np.square(X - centerfilter) / (2 * np.square(sigma_gfilter)))

    return gfilter, sigma_gfilter


'''
########################################################################################
#Windows LAYOUTS
Building frames for main windows
########################################################################################
'''

# LAYOUT INTERFEROMETER IMAGE
layout_frame_ImgSample = [
    [sg.Image(size=size, background_color='black', key='image1', expand_x=True)],
    [sg.Input(expand_x=True, disabled=True, key='file1', visible='True')],
    [sg.Button('Open File(s)', font='Arial 10 bold'),
     sg.Button('Rotate (°)', visible=True, font='Arial 10 bold', disabled=True),
     sg.Input('0', size=(5, 1), key='-DEGREE-', enable_events=True),
     sg.Text('Original Size (w,h):'),
     sg.Text('', key='-scale1-')],
]
# LAYOUT REFERENCE IMAGE
layout_frame_ImgReference = [
    [sg.Image(size=size2, background_color='black',
              key='image2', enable_events=True)],
    [sg.Input(expand_x=True, disabled=True, key='file2', visible='True')],
    [sg.Button('Open Ref.', font='Arial 10 bold')],
]

# LAYOUT SELECT AREA OPTIONS
layout_area_selection = [
    [sg.Text('X Coord')],
    [sg.Spin([i for i in range(minvalue_x, maxvalue_x + 1)], initial_value=begin_x, key='-BEGIN_X-', size=(5, 1),
             enable_events=True),
     sg.Spin([i for i in range(minvalue_x, maxvalue_x + 1)], initial_value=end_x, key='-END_X-', size=(5, 1),
             enable_events=True)],

    [sg.Text('Y Coord')],
    [sg.Spin([i for i in range(minvalue_y, maxvalue_y + 1)], initial_value=begin_y, key='-BEGIN_Y-', size=(5, 1),
             enable_events=True),
     sg.Spin([i for i in range(minvalue_y, maxvalue_y + 1)], initial_value=end_y, key='-END_Y-', size=(5, 1),
             enable_events=True)],
    [sg.Text('BG Phase Fit. (%): ')],
    [sg.Input(base_ref, size=(6, 1), key='-base_ref-', enable_events=True)],
]

# LAYOUT INPUT GAS AND RADIATION PARAMETERS
layout_input_parameters = [
    [sg.Text('Scaling Factor (µm/pixel):         '),
     sg.Input(factor, size=(6, 1), key='-factor-', enable_events=True)],
    [sg.Text('Unc. Scaling Factor (µm/pixel): '),
     sg.Input(unc_factor, size=(6, 1), key='-unc_factor-', enable_events=True)],
    [sg.Text('Fringes Orientation:      '),
     sg.Combo(['vertical', 'horizontal'], default_value='horizontal', key='-combofringe-')],
    [sg.Text('Laser Wavelength λ (nm):          '),
     sg.Input(lambda0, size=(6, 1), key='-lambda0-', enable_events=True)],
    [sg.Text('Wavelength FWHM Δλ (nm):     '),
     sg.Input(unc_lambda0, size=(6, 1), key='-unclambda0-', enable_events=True)],
    [sg.Text('Gas:'),
     sg.Combo(['H2', 'N2', 'He', 'Ar', '--'], default_value='N2', key='-combogas-', enable_events=True),
     sg.Text(' Polarizab. (Å³): '),
     sg.Input(polargas, size=(6, 1), key='-polargas-', enable_events=True)],
]
# LAYOUT INPUT MEASUREMENT PARAMETERS
layout_analysis_parameters = [

    [sg.Text('Filter Frequency   ν\u2080 (pixel):      '),
     sg.Spin([i for i in range(minvalue_x, maxvalue_x + 1)], initial_value=0, size=(5, 1), key='-centerfilter-',
             enable_events=True)],
    [sg.Text('Filter Range  Δν (pixel):             '),
     sg.Spin([i for i in range(0, maxvalue_x // 2)], initial_value=0, size=(5, 1), key='-sigma_gfilter-',
             enable_events=True)],
    [sg.Text('Gaussian Blur σ (pixel):            '),
     sg.Spin([i for i in range(0, maxvalue_x // 2)], initial_value=sigma_gblur, size=(5, 1), key='-sigma_gblur-',
             enable_events=True)],
    [sg.Text('')],
    [sg.Text('Axisymmetric Orientation:'),
     sg.Combo(['vertical', 'horizontal','none'], default_value='vertical', enable_events = True, key='-comboaxisymm-')],
    [sg.Text('Axisymetric Position (pixel):     ', key = '-text_axissymm-'),
     sg.Spin([i for i in range(minvalue_x, maxvalue_x + 1)], initial_value=0, size=(5, 1), key='-axisymm_pos-',
             enable_events=True)],
]
# LAYOUT FRAME OF ALL INPUT OPTIONS
layout_frame_Options = [
    [sg.Frame('Select Area', layout_area_selection, size=(130, 210), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 10 bold'),
     sg.Frame('Input Parameters', layout_input_parameters, size=(260, 210), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 10 bold'),
     sg.Frame('Analysis Parameters', layout_analysis_parameters, size=(260, 210), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 10 bold')],
]
# LAYOUT FRAME LEFT - INPUTS
layout_frame_ImagesL = [
    [sg.Frame("Interferogram (Target)", layout_frame_ImgSample, size=(440, 430), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 12 bold')],
]
# LAYOUT FRAME RIGHT - INPUTS AND APPLY
layout_frame_ImagesR = [
    [sg.Frame("Interferogram (Ref.)", layout_frame_ImgReference, size=(258, 258), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 12 bold')],
    [sg.Button('Analyse Data', size=(30, 6), font='Arial 12 bold', disabled=True, button_color='black')],
    [sg.Button('Clear', size=(30, 2), button_color='gray', font='Arial 10 bold')]
]
# lAYOUT GLOBAL INPUTS
layout_frame_Images = [
    [sg.Column(layout_frame_ImagesL, element_justification='c'),
     sg.Column(layout_frame_ImagesR, element_justification='c')],
    [sg.Frame("Options", layout_frame_Options, size=(698, 210), title_location=sg.TITLE_LOCATION_TOP_LEFT,
              font='Arial 12 bold')],
]
# LAYOUT 1D OPTIONS PLOT
layout_frame_plot1D = [
    [sg.Checkbox('Prof. 1 (µm)', default=False, key='-checkpos1-'),
     sg.Input(pos1, size=(4, 1), key='-pos1-', enable_events=True),
     sg.Checkbox('Prof. 2 (µm)', default=False, key='-checkpos2-'),
     sg.Input(pos2, size=(4, 1), key='-pos2-', enable_events=True),
     sg.Checkbox('Prof. 3 (µm)', default=False, key='-checkpos3-'),
     sg.Input(pos3, size=(4, 1), key='-pos3-', enable_events=True)],
    [sg.Slider(key='sliderh', range=(0, 490), orientation='h', size=(400, 20), default_value=0,
               enable_events=True)]
]
# LAYOUT STAGES OF THE TREATMENT
layout_frame_Steps = [
    [sg.Radio('Frequency \nDomain', "RADIO1", default=False, key='fftradio'),
     sg.Radio('Gaussian \nFilter', 'RADIO1', default=False, key='filterradio'),
     sg.Radio('Acc. \nPhase-shift', "RADIO1", default=True, key='phaseradio'),
     sg.Radio('Radial \nPhase-shift', "RADIO1", default=False, key='abelradio', enable_events=True),
     sg.Radio('Density \nProfile', "RADIO1", default=False, key='densradio')],
]
# LAYOUT GLOBAL OUTPUTS
layout_frame_Result = [
    [sg.Frame('Stages', layout_frame_Steps, size=(490, 70), title_location=sg.TITLE_LOCATION_TOP_LEFT,
              key='framesteps', font='Arial 10 bold')],
    [sg.Button('1D Profile', size=(16, 1), font='Arial 10 bold', disabled=True),
     sg.Button('2D Profile', size=(16, 1), font='Arial 10 bold', disabled=True),
     sg.Checkbox('Unc. Measurement', default=False, key='-checkstd-'),
     ],
    [sg.Canvas(key='canvasabel', size=(490, 400), background_color='black')],
    [sg.Button('Save Plot', size=(16, 2), disabled=True, font='Arial 10 bold'),
     sg.Button('Save Data', size=(16, 2), disabled=True, font='Arial 10 bold'),
     sg.Text('Colormap :'),
     sg.Combo(['default', 'rainbow', 'plasma', 'viridis'], default_value='default', key='-cmapcombo-', enable_events=True)
     ],
    [sg.Frame('1D Profile (Axisymmetry)', layout_frame_plot1D, size=(490, 100),
              title_location=sg.TITLE_LOCATION_TOP_LEFT,
              visible=False, key='frame1d', font='Arial 10')],

]
# lAYOUT GLOBAL
layout = [
    [sg.Frame("Interferograms", layout_frame_Images, size=(700, 700), font='Arial 12 bold'),
     sg.Frame("Target Profile", layout_frame_Result, size=(500, 700), title_location=sg.TITLE_LOCATION_TOP,
              font='Arial 12 bold')],
]
######################################################################################################################
window = sg.Window("Interferogram Analysis - GJ (Version 2.1.1)", layout, margins=(1, 1), finalize=True)
######################################################################################################################
'''
####################################################################################################
#WINDOWS EVENT
####################################################################################################
'''
mouse_draw = False
mclick = False
window['image1'].bind('<ButtonPress-1>', '-click1-')
window['image1'].bind('<ButtonRelease-1>', '-click2-')
window['-BEGIN_X-'].bind('<FocusOut>', 'FocusOut')
window['-END_X-'].bind('<FocusOut>', 'FocusOut')
window['-BEGIN_Y-'].bind('<FocusOut>', 'FocusOut')
window['-END_Y-'].bind('<FocusOut>', 'FocusOut')

while True:
    event, values = window.read()
    # Removing temp files when the main window is closed
    if event == sg.WINDOW_CLOSED:
        if '_temp.png' in path1:
            os.remove(path1)
            os.remove(path2)
        break

    if event == 'Clear':
        if '_temp.png' in path1:
            os.remove(path1)
            os.remove(path2)
        window['Rotate (°)'].update(disabled=True)
        window['-DEGREE-'].update(visible=True)
        window['Analyse Data'].update(disabled=True)

        # Disable specific buttons and frames for 2D analysis
        window['Save Data'].update(disabled=True)
        window['Save Plot'].update(disabled=True)
        window['frame1d'].update(visible=False)
        window['2D Profile'].update(disabled=True)
        window['1D Profile'].update(disabled=True)
        window['Rotate (°)'].update(disabled=True)
        window['-DEGREE-'].update(visible=True)
        window['Analyse Data'].update(disabled=True)

        # Reset values
        path1 = ''
        path2 = ''
        window['image1'].update(size=size, data='')
        window['image2'].update(size=size2, data='')
        window['file1'].update(path1)
        window['file2'].update(path1)
        window['-centerfilter-'].update('0')
        window['-sigma_gfilter-'].update('0')
        window['-axisymm_pos-'].update('0')
        # Cleaning plots
        try:
            fig_canvas_agg.get_tk_widget().forget()
        except:
            plt.close('all')

    #######################################################################
    if event == 'image1-click1-':
        e = window['image1'].user_bind_event
        if mclick == False:
            window["-BEGIN_X-"].update(f'{e.x}')
            window["-BEGIN_Y-"].update(f'{e.y}')
            centerfilter = 0
            axis_pos = 0
            window['-centerfilter-'].update('0')
            window['-sigma_gfilter-'].update('0')
            window['-axisymm_pos-'].update('0')
            mclick = True
        else:
            window["-END_X-"].update(f'{e.x}')
            window["-END_Y-"].update(f'{e.y}')
            mclick = False

    if event == 'image1-click2-' and mclick == False:
        if path1 != '' and path2 != '':
            apply_drawing(values, window)
    ########################################################################
    # SYMMETRIC OR ASSYMETRIC TARGET ANALYSIS
    if event == '-comboaxisymm-':
        if values['-comboaxisymm-'] == 'none':
            window['-text_axissymm-'].update('Target Thickness (µm):            ')
            window['abelradio'].update(disabled=True)
        else:
            window['-text_axissymm-'].update('Axisymetric Position (pixel):     ')
            window['abelradio'].update(disabled=False)
    ########################################################################
    # OPEN INTERFEROGRAM IMAGE
    elif event == "Open File(s)":
        if '_temp.png' in path1:
            try:
                os.remove(path1)
            except UnidentifiedImageError:
                continue

        path_files = sg.popup_get_file("", no_window=True, multiple_files=True)
        if path_files:
            path1 = path_files[0]
        else:
            path1 = ''
        window['file1'].update(path1)
        # No file open
        if path1 == '':
            continue
        # create PNG files from files SNP
        # Note: files with SNP extension must be converted to PNG for algorithm analysis
        if '.snp' in path1:
            path_files_snp = path_files
            path_files = []
            for ipath in path_files_snp:
                databinary = getBinaryData(ipath)
                data0 = np.flip(databinary[60:60 + 720 * 576])
                originaltgtsnp = Image.new(mode='L', size=(720, 576))
                originaltgtsnp.putdata(data0)

                ipath = ipath.replace('.snp', '_temp.png')

                originaltgtpng = originaltgtsnp.save(ipath)
                if len(path_files) == 0:
                    path1 = ipath
                    window['file1'].update(path1)

                path_files.append(ipath)

        try:
            # Open Files
            originaltgt0 = []
            for i in range(0, len(path_files)):
                originaltgt0.append(Image.open(path_files[i]))
            apply_drawing(values, window)

        except UnidentifiedImageError:
            continue

        # scale 1: scale for interferogram image
        w, h = originaltgt0[0].size
        scale = (width / w), (height / h)

        if scale != 1:
            im1 = originaltgt0[0].resize(size)
        else:
            im1 = originaltgt0[0]
        data1 = image_to_data(im1)

        window['image1'].update(data=data1, size=size)
        window['-scale1-'].update(originaltgt0[0].size)
        centerfilter = 0
        axis_pos = 0
        window['-centerfilter-'].update('0')
        window['-axisymm_pos-'].update('0')
        # Enable buttons
        if path1 != '' and path2 != '':
            window['Rotate (°)'].update(disabled=False)
            window['-DEGREE-'].update(visible=True)
            window['Analyse Data'].update(disabled=False)

    ########################################################################
    # OPEN REFERENCE FILE
    elif event == "Open Ref.":
        if '_temp.png' in path2:
            try:
                os.remove(path2)
            except UnidentifiedImageError:
                continue

        path_files2 = sg.popup_get_file("", no_window=True, multiple_files=True)
        if path_files2:
            path2 = path_files2[0]
        else:
            path2 = ''
        window['file2'].update(path2)
        # No file
        if path2 == '':
            continue
        # Create PNG files from files SNP
        if '.snp' in path2:
            path_files_snp2 = path_files2
            path_files2 = []
            for ipath2 in path_files_snp2:
                databinary2 = getBinaryData(ipath2)
                data02 = np.flip(databinary2[60:60 + 720 * 576])
                originalrefsnp = Image.new(mode='L', size=(720, 576))
                originalrefsnp.putdata(data02)

                ipath2 = ipath2.replace('.snp', '_temp.png')

                originalrefpng = originalrefsnp.save(ipath2)
                if len(path_files2) == 0:
                    path2 = ipath2
                    window['file2'].update(path2)

                path_files2.append(ipath2)
        try:
            # Open Files
            originalref0 = []
            for i in range(0, len(path_files2)):
                originalref0.append(Image.open(path_files2[i]))
            apply_drawing(values, window)

        except UnidentifiedImageError:
            continue
        w2, h2 = originalref0[0].size
        scale2 = (width2 / w2), (height2 / h2)

        if scale2 != 1:
            im2 = originalref0[0].resize(size2)
        else:
            im2 = originalref0[0]

        data2 = image_to_data(im2)
        window['image2'].update(data=data2, size=size2)
        centerfilter = 0
        axis_pos = 0
        window['-centerfilter-'].update('0')
        window['-axisymm_pos-'].update('0')
        # Enable buttons
        if path1 != '' and path2 != '':
            window['Rotate (°)'].update(disabled=False)
            window['-DEGREE-'].update(visible=True)
            window['Analyse Data'].update(disabled=False)

    ########################################################################
    # BUTTON COORD AREA
    elif event == '-BEGIN_X-' or event == '-END_X-':
        apply_drawing(values, window)
        centerfilter = 0
        axis_pos = 0
        window['-centerfilter-'].update('0')
        window['-axisymm_pos-'].update('0')
    elif event == '-BEGIN_X-FocusOut' or event == '-END_X-FocusOut':
        centerfilter = 0
        axis_pos = 0
        window['-centerfilter-'].update('0')
        window['-axisymm_pos-'].update('0')
    elif event == '-BEGIN_Y-' or event == '-END_Y-':
        apply_drawing(values, window)
        centerfilter = 0
        axis_pos = 0
        window['-centerfilter-'].update('0')
        window['-axisymm_pos-'].update('0')
    elif event == '-BEGIN_Y-FocusOut' or event == '-END_Y-FocusOut':
        centerfilter = 0
        axis_pos = 0
        window['-centerfilter-'].update('0')
        window['-axisymm_pos-'].update('0')
    #########################################################################
    # BUTTON ROTATE
    elif event == 'Rotate (°)':
        apply_drawing(values, window)
        centerfilter = 0
        axis_pos = 0
        window['-centerfilter-'].update('0')
        window['-axisymm_pos-'].update('0')
    '''
    #######################################################################
    # BUTTON APPLY - main event of window
    In this event will be apply the treatment of interferogram image to generate
    the data of gas profile.  
    #######################################################################
    '''
    if event == 'Analyse Data' or event == '-axisymm_pos-' or event == '-sigma_gblur-' \
            or event == '-centerfilter-' or event == '-sigma_gfilter-':
        apply_drawing(values, window)
        # Cleaning plots
        try:
            fig_canvas_agg.get_tk_widget().forget()
        except:
            plt.close('all')
        # Input datas
        h_prof = -1.0
        try:
            # get rectangle coord.
            begin_x = int(get_value("-BEGIN_X-", values) / scale[0])
            begin_y = int(get_value("-BEGIN_Y-", values) / scale[1])
            end_x = int(get_value("-END_X-", values) / scale[0])
            end_y = int(get_value("-END_Y-", values) / scale[1])

            if begin_x > end_x:
                begin_x = int(get_value("-END_X-", values) / scale[0])
                end_x = int(get_value("-BEGIN_X-", values) / scale[0])
                window["-END_X-"].update(int(begin_x * scale[0]))
                window["-BEGIN_X-"].update(int(end_x * scale[0]))

            elif begin_x == end_x:
                end_x = begin_x + 100
                window["-END_X-"].update(value=int(end_x * scale[0]))
                window["-BEGIN_X-"].update(value=int(begin_x * scale[0]))

            if begin_y > end_y:
                begin_y = int(get_value("-END_Y-", values) / scale[1])
                end_y = int(get_value("-BEGIN_Y-", values) / scale[1])
                window["-END_Y-"].update(value=int(begin_y * scale[1]))
                window["-BEGIN_Y-"].update(value=int(end_y * scale[1]))

            elif begin_y == end_y:
                end_y = begin_y + 100
                window["-END_Y-"].update(value=int(end_y * scale[1]))
                window["-BEGIN_Y-"].update(value=int(begin_y * scale[1]))

            # get angle to image rotation
            rotate_degree = -float(get_value('-DEGREE-', values))
            # get base ref %
            base_ref = float(get_value('-base_ref-', values)) / 100
            # get conversion factor in meters/pixel
            factor = float(get_value('-factor-', values)) * 1e-6
            unc_factor = float(get_value('-unc_factor-', values)) * 1e-6
            # Manual definition of the filter position
            centerfilter = int(get_value("-centerfilter-", values))
            sigma_gfilter = int(get_value('-sigma_gfilter-', values))
            # sigma value of gaussian blur
            sigma_gblur = int(get_value('-sigma_gblur-', values))
            # Manual axisymm position
            axis_pos = int(get_value('-axisymm_pos-', values)) # pixel
            # heat gas constant
            alpha_gas = float(get_value('-polargas-', values)) * 1e-30
            # Wavelength laser
            lambda0 = float(get_value('-lambda0-', values)) * 1e-9  # in meters
            unc_lambda0 = float(get_value('-unclambda0-', values)) * 0.588705 * 1e-9  # sigma in meters
            # Manual axisymm position
            if values['-comboaxisymm-'] == 'none':
                axis_pos = 0
                thick_tgt =  float(get_value('-axisymm_pos-', values))*1e-6/factor
                std_thick_tgt = thick_tgt/factor*unc_factor
            else:
                axis_pos = int(get_value('-axisymm_pos-', values))  # pixel
                thick_tgt = 0

        except:
            sg.popup_error(f"WARNING: Data fields must have numerical values! ")
            continue

        newcmp = cm.get_cmap(values['-cmapcombo-'], 256)

        # Rotate Image
        if rotate_degree != 0:
            originalref = rotate(originalref0[0], rotate_degree, reshape = False)
        else:
            originalref = originalref0[0]
        # Input ref. array from ref. image
        array_ref = np.asarray(originalref)

        # Input original files
        phasemaps = []
        tgt_dens, tgt_abelphasemap, tgt_phasemap = [], [], []
        std_phasemap, std_abelmap, std_tgt_dens = [], [], []

        for j in range(0, len(path_files)):
            if len(path_files) == len(path_files2):
                if rotate_degree != 0:
                    originaltgt = rotate(originaltgt0[j], rotate_degree, reshape=False)
                    originalref = rotate(originalref0[j], rotate_degree, reshape=False)
                else:
                    originalref = originalref0[j]
                    originaltgt = originaltgt0[j]
            else:
                if rotate_degree != 0:
                    originaltgt = rotate(originaltgt0[j], rotate_degree, reshape=False)
                else:
                    originaltgt = originaltgt0[j]

            array_tgt = np.asarray(originaltgt)
            array_ref = np.asarray(originalref)

            if np.ndim(array_tgt) == 3:
                # Slice image with 3 channels:only one channel is used to interferogram treatment
                intref0 = array_ref[:, :, 0]
                inttgt0 = array_tgt[:, :, 0]
            else:
                intref0 = array_ref[:, :]
                inttgt0 = array_tgt[:, :]

            # select area of image
            intref = intref0[begin_y:end_y, begin_x:end_x]
            inttgt = inttgt0[begin_y:end_y, begin_x:end_x]


            try:
                # Apply Fast Fourier Transform on interferogram data arrays
                fftref = np.fft.fft2(intref)  # ref. interferogram
                ffttgt = np.fft.fft2(inttgt)  # gas interferogram
                # Defining line or row to apply gaussian filter
                fftmap = np.log(np.abs(ffttgt))
                nlmap, nrmap = np.shape(fftmap)
            except:
                fftref = np.zeros(np.shape(intref))  # ref. interferogram
                ffttgt = np.zeros(np.shape(intref))  # gas interferogram
                # Defining line or row to apply gaussian filter
                fftmap = np.zeros(np.shape(intref))
                nlmap, nrmap = np.shape(fftmap)

            rotate_degree = 0
            '''
            # Authomatic definition of the gaussian filter position: 
            this position are defined like the line or column (Vertical or horizontal fringes) with more intensity pixel
            values. This way, this positions are defined using the maximum value of horizontal/vertical pixels sum, 
            depending on fringes orientation.
            Note: Case this filter position is not found, the process is interrupted and the user must select another
            file or another area of interferogram.
            '''
            fftmap = (fftmap - np.min(fftmap)) * np.ones(np.shape(fftmap)) / (np.max(fftmap) - np.min(fftmap))

            if centerfilter == 0:
                summap = []
                # Fringe Orientation: HORIZONTAL or VERTICAL
                # sum of array rows (vertical)
                if values['-combofringe-'] == 'vertical':
                    summap = np.sum(fftmap, axis=0) - als(np.sum(fftmap, axis=0))
                else:
                    # sum of array lines (horizontal)
                    summap = np.sum(fftmap, axis=1) - als(np.sum(fftmap, axis=1))

                try:
                    # Defining point of gaussian filter application using max value of horizontal/vertical pixels sum
                    filterpoints, _ = find_peaks(summap, height=0.5 * np.max(summap), width = 1)
                    # Range of gaussian filter
                    filterspoints_widths = (peak_widths(summap, filterpoints, rel_height=0.5)[0])

                except:
                    filterpoints = []
                    filterspoints = []

                if len(filterpoints) == 0:
                    sg.popup(f"WARNING: Unable to apply the Fast Fourier Transform to the selected image!")
                    continue
                if values['-combofringe-'] == 'horizontal':
                    # filter range is equal to FWHM of signal of summaps
                    if filterpoints[0] <= 5:
                        centerfilter = filterpoints[1]
                        f_range = 2 * int(filterspoints_widths[1])
                    else:
                        centerfilter = filterpoints[0]
                        f_range = 2 * int(filterspoints_widths[0])

                elif values['-combofringe-'] == 'vertical':

                    if filterpoints[len(filterpoints) - 1] >= nrmap - 5:
                        centerfilter = filterpoints[len(filterpoints) - 2]
                        f_range = 2 * int(filterspoints_widths[len(filterpoints) - 2])
                    else:
                        centerfilter = filterpoints[len(filterpoints) - 1]
                        f_range = 2 * int(filterspoints_widths[len(filterpoints) - 1])

            window['-centerfilter-'].update(str(centerfilter))

            # Creating filter array from null array
            gfilter, sigma_gfilter = func_gfilter(fftref, centerfilter, f_range, sigma_gfilter)

            # Applying Inverse FFT in resultant array obtained after use of the gaussian filter on FFT arrays
            ifftref = np.fft.ifft2(gfilter * fftref)
            iffttgt = np.fft.ifft2(gfilter * ffttgt)

            # Creating Phase Maps arrays by subtracting the arguments of IFFT arrays
            phasemaps = (np.angle(iffttgt) - np.angle(ifftref))
            # Unwrap phase:
            uwphasemap = unwrap_phase(phasemaps)

            '''##########################################################################################
            DEFINING STANDARD DEVIATION:
            The standard deviation is calculated from fringes intensity distribution, fringes widths and 
            fringes displacement.
            '''  #########################################################################################
            # Estimating displacement (vertical and horizontal) between interferograms
            try:
                disp_xy, _, _ = phase_cross_correlation(inttgt0, intref0, upsample_factor=100)
            except:
                disp_xy = np.array([0.0,0.0])

            if values['-combofringe-'] == 'vertical':
                disp = np.absolute(disp_xy[1])
                # Creating 2D array for fringes width distribution
                dist_fw = fringes_width(intref)

            if values['-combofringe-'] == 'horizontal':
                disp = np.absolute(disp_xy[0])
                # Creating 2D array for fringes width distribution
                dist_fw = np.transpose(fringes_width(np.transpose(intref)))

            '''
            NOTE: During our algorithm tests we verify some computational artefacts. 
            These artifacts are detected only in the multiplication of the intensity distributions. 
            To correct this error we add a baseline line over data. The baseline has a value equal 
            to 0.1% of the lesser intensity.  
            '''
            distI1 = np.asfarray(intref) - np.min(intref)
            distI2 = np.asfarray(inttgt) - np.min(inttgt)
            func_distI = np.sqrt((np.mean(distI1) * (distI1 + distI2))/(2*distI1*distI2))

            try:
                std_phasemap_i = ((np.pi * disp) / (2 * dist_fw)) * (np.nan_to_num(func_distI, posinf = np.median(func_distI)))

            except:
                std_phasemap_i = np.zeros(np.shape(intref))
            '''
            ################################################################################
            Applying Inverse Abel Transform (IAT):
            The IAT is applied using Dash Onion Peeling algorithm from PyAbel. To apply its library correctly is necessary
            to define a axis symmetric in image (Horizontal or Vertical) and it is defined from more intensity pixel range. 
            So, the image is cut according axissymetric.
            NOTE: the Abel transform is always performed around the vertical axis, so when the image have horizontal
            axissymmetry the matrix must be transposed.

            '''

            # Transpose Matrix for Horizontal Axissmetry
            if values['-comboaxisymm-'] == 'horizontal':
                uwphasemap = np.transpose(uwphasemap)
                std_phasemap_i = np.transpose(std_phasemap_i)

            # Remove rising background of PIL
            nlines, nrows = np.shape(uwphasemap)

            baselinemap, std_blmap = baseline2D(uwphasemap, base_ref)

            resultphase = (uwphasemap - baselinemap) - np.min(uwphasemap - baselinemap) * np.ones(np.shape(uwphasemap))

            '''
            Breakpoint()
            #.dat Files for manuscript
            uwphasemap = (gaussian_filter(uwphasemap, sigma=sigma))
            baselinemap = (gaussian_filter(baselinemap, sigma=sigma))
            resultphase = (gaussian_filter(resultphase, sigma=sigma))

            np.savetxt('phasemap+bg pos00.dat', uwphasemap, fmt='%.3e')
            np.savetxt('bg pos00.dat', baselinemap, fmt='%.3e')
            np.savetxt('phasemap-bg pos00.dat', resultphase, fmt='%.3e')
            '''

            # First contribution: Experimental + baseline
            std_phasemap_i = np.sqrt(np.square(std_blmap * np.ones(np.shape(uwphasemap))) + np.square(std_phasemap_i))

            # Apply gaussian blur
            phasemap_corr = (gaussian_filter(resultphase, sigma=sigma_gblur))
            # Std without gaussian blur
            std_resultphase = ((np.ones(np.shape(resultphase))).T * np.std(resultphase, axis=1)).T

            # Std caused by Gaussian blur
            if sigma_gblur == 0:
                std_phasemap_corr = np.zeros(np.shape(resultphase))
            else:
                std_phasemap_corr = np.abs(gaussian_filter(resultphase, sigma=1) - phasemap_corr)

            std_phasemap_corr = np.sqrt(np.square(std_phasemap_corr) + np.square(std_phasemap_i) +
                                        np.square(std_resultphase))

            # Data processing for asymmetrical targets
            if values['-comboaxisymm-'] == 'none':

                phasemap_symm = np.zeros(np.shape(uwphasemap))
                std_phasemap_symm = np.zeros(np.shape(uwphasemap))
                std_abelmap_2 = np.zeros(np.shape(uwphasemap))
                phase_abel = np.zeros(np.shape(uwphasemap))
                norm_phasemap = np.zeros(np.shape(uwphasemap))

                if thick_tgt == 0:
                    thick_tgt = 1
                    window['-axisymm_pos-'].update(thick_tgt)

                '''
                ########################################################################################
                Calculating refraction index and plasma electronic density from IAT phasemap.        
                '''
                # Calculating index refraction from IAT of phasemap
                n_index0 = (1 + (phasemap_corr * lambda0) / (2 * np.pi * thick_tgt*factor))
                # Cutting border of images due the computational artefacts generated by IAT and problems with no symmetric images
                #n_index = n_index0[:, side_cut: -side_cut]
                # Calculating gas density
                const_gas = 3 / (4 * np.pi * alpha_gas)

                try:
                    tgt_dens_i = const_gas * ((n_index0 ** 2 - 1) / (n_index0 ** 2 + 2)) * 1e-6  # cm-3

                except:
                    tgt_dens_i = np.zeros(np.shape(uwphasemap))

                # new matrix size for plot
                rangeh, rangev = np.shape(uwphasemap)

                '''
                CALCULATION TOTAL STANDARD DEVIATION FROM:
                1. Measurement of interferogram
                2. Inverse Abel Transform
                3. FWHM Laser wavelength
                '''
                dN_phase = const_gas * (lambda0 / (2 * np.pi * thick_tgt * factor)) * np.ones(np.shape(uwphasemap))  # rad

                # Contribution 3: laser wavelength
                dN_lambda = const_gas * (phasemap_corr / (2 * np.pi * thick_tgt * factor))  # rad/m

                # Contribution 4: Scaling factor unc.
                dN_factor = const_gas * (phase_abel * lambda0 / (2 * np.pi * factor ** 2))

                std_tgt_dens_i = np.sqrt(np.square(dN_phase * std_phasemap_corr) + np.square(dN_lambda * unc_lambda0) +
                                         np.square(dN_factor * unc_factor)) * 1e-6

            #####################################################################################################
            # Data processing for symmetrical targets
            else:

                # Apply gaussian filter to define the region with more intensity pixel value
                # Define region with more intensity pixel - position x and y
                if axis_pos == 0:
                    try:
                        cline, crow = np.where(np.abs(phasemap_corr) >= np.max(np.abs(phasemap_corr)) * 0.98)
                        cx = int(np.median(crow))
                    except:
                        cx = np.unravel_index(np.argmax(np.abs(phasemap_corr), axis=None), phasemap_corr.shape)[1]

                if axis_pos != 0:
                    if np.shape(phasemaps) == np.shape(uwphasemap):
                        cx = axis_pos
                    else:
                        cx = nrows - axis_pos

                # If the region not found, set symmetric point like half image
                if math.isnan(cx) == True:
                    cx = int(nrows / 2)

                if cx <= int(nrows / 2):
                    vert_lim = int(2 * cx + 1)
                    phasemap_symm = phasemap_corr[:, 0:vert_lim]
                    std_phasemap_symm = std_phasemap_corr[:, 0:vert_lim]

                # If right-side of image is more width
                if cx > int(nrows / 2):
                    vert_lim = int(2 * cx - nrows)
                    phasemap_symm = phasemap_corr[:, vert_lim:]
                    std_phasemap_symm = std_phasemap_corr[:, vert_lim:]

                # Creating axisymetric line in plot
                if np.shape(phasemaps) == np.shape(uwphasemap):
                    axis_pos = cx
                else:
                    axis_pos = nrows - cx

                window['-axisymm_pos-'].update(axis_pos)

                try:

                    raxis = np.abs(np.arange(-np.shape(phasemap_symm)[1]/ 2, np.shape(phasemap_symm)[1]/ 2, 1))

                    # Applying inverse Abel Transform
                    phase_abel0 = abel.Transform((phasemap_symm), symmetry_axis=0, direction='inverse',
                                                 method='onion_peeling').transform

                    std_abel0 = std_phasemap_symm / (np.pi * np.max(raxis) * np.sqrt(1 - raxis / np.max(raxis)))

                except:

                    phase_abel0 = np.zeros(np.shape(phasemap_corr[:, 0:vert_lim]))
                    std_abel0 = np.zeros(np.shape(phasemap_corr[:, 0:vert_lim]))
                    sg.popup_error(f"WARNING: Unable to apply the Abel transform to the selected image! ")
                    continue

                '''
                ############################################################################################
                Calculating std from Abel Transform:
                The std is calculated using deviation of mormalized phasemap and normalized IAT phasemap  
                '''
                rangeh0, rangev0 = np.shape(phase_abel0)
                # side cut to Abel (%)
                side_cut = int(0.05 * rangev0)

                phase_abel = phase_abel0[:, side_cut: -side_cut]
                norm_phasemap = np.zeros(np.shape(phase_abel))
                for k in range(0, rangeh0):
                    try:
                        norm_phasemap[k] = phasemap_symm[k, side_cut: -side_cut] \
                                           * np.max(abs(phase_abel[k, :])) \
                                           / np.max(abs(phasemap_symm[k, side_cut: -side_cut]))
                    except:
                        norm_phasemap[k] = np.zeros(np.shape(phase_abel[k, :]))

                std_abelmap_2 = np.sqrt(np.square(phase_abel - norm_phasemap))

                '''
                ########################################################################################
                Calculating refraction index and plasma electronic density from IAT phasemap.        
                '''
                # Calculating index refraction from IAT of phasemap
                n_index0 = (1 + (phase_abel0 * lambda0) / (2 * np.pi * factor))
                # Cutting border of images due the computational artefacts generated by IAT and problems with no symmetric images
                n_index = n_index0[:, side_cut: -side_cut]
                # Calculating gas density
                const_gas = 3 / (4 * np.pi * alpha_gas)

                try:
                    tgt_dens_i = const_gas * ((n_index ** 2 - 1) / (n_index ** 2 + 2)) * 1e-6  # cm-3

                except:
                    tgt_dens_i = np.zeros(np.shape(phase_abel))

                # new matrix size for plot
                rangeh, rangev = np.shape(phase_abel)

                '''
                CALCULATION TOTAL STANDARD DEVIATION FROM:
                1. Measurement of interferogram
                2. Inverse Abel Transform
                3. FWHM Laser wavelength
                '''
                # Contribution 1: measurement interferogram
                std_phase1 = abs(std_abel0[:, side_cut: -side_cut] / factor)  # rad/um

                # Contribution 2: Phase
                std_phase2 = abs(std_abelmap_2 / factor)  # rad/m

                std_phase = np.sqrt(np.square(std_phase1) + np.square(std_phase2))  # rad/m

                dN_phase = const_gas * (lambda0 / (2 * np.pi)) * np.ones(np.shape(phase_abel))  # rad

                # Contribution 3: laser wavelength
                dN_lambda = const_gas * (phase_abel / (2 * np.pi * factor))  # rad/m

                # Contribution 4: Scaling factor unc.
                dN_factor = const_gas * (phase_abel * lambda0 / (2 * np.pi * factor ** 2))

                std_tgt_dens_i = np.sqrt(np.square(dN_phase * std_phase) + np.square(dN_lambda * unc_lambda0) +
                                         np.square(dN_factor * unc_factor)) * 1e-6

            if values['-comboaxisymm-'] == 'horizontal':
                phasemap_corr = np.transpose(phasemap_corr)
                phase_abel = np.transpose(phase_abel)
                std_abelmap_2 = np.transpose(std_abelmap_2)
                std_phase = np.transpose(std_phase)
                std_phase1 = np.transpose(std_phase1)
                std_tgt_dens_i = np.transpose(std_tgt_dens_i)
                tgt_dens_i = np.transpose(tgt_dens_i)
                std_phasemap_corr = np.transpose(std_phasemap_corr)
                norm_phasemap = np.transpose(norm_phasemap)

            tgt_phasemap.append(phasemap_corr)
            std_phasemap.append(std_phasemap_corr)

            tgt_abelphasemap.append(phase_abel)
            std_abelmap.append(std_abelmap_2)

            tgt_dens.append(tgt_dens_i)
            std_tgt_dens.append(std_tgt_dens_i)

        # BUILDING 2D ARRAYS RESULTS FOR:
        if len(tgt_phasemap) > 1:  # Many files
            # PHASEMAP
            tgt_phasemap_mean = mean_maps(tgt_phasemap)
            std_phasemap_mean = np.sqrt(np.square(mean_maps(std_phasemap)) + \
                                        np.square(std_maps(tgt_phasemap, tgt_phasemap_mean)))
            # INV. ABEL TRANSF. MAP
            tgt_abelmap_mean = mean_maps(tgt_abelphasemap)
            std_abelmap_mean = np.sqrt(np.square(mean_maps(std_abelmap)) + \
                                       np.square(std_maps(tgt_abelphasemap, tgt_abelmap_mean)))
            # PLASMA DENSITY
            tgt_dens_mean = mean_maps(tgt_dens)
            std_dens_mean = np.sqrt(np.square(mean_maps(std_tgt_dens)) + \
                                    np.square(std_maps(tgt_dens, tgt_dens_mean)))
        else:
            try:
                # PHASEMAP
                tgt_phasemap_mean = tgt_phasemap[0]
                std_phasemap_mean = std_phasemap[0]
                # INV. ABEL TRANSF. MAP
                tgt_abelmap_mean = tgt_abelphasemap[0]
                std_abelmap_mean = std_abelmap[0]
                # PLASMA DENSITY
                tgt_dens_mean = tgt_dens[0]
                std_dens_mean = std_tgt_dens[0]
            except:
                # PHASEMAP
                tgt_phasemap_mean = np.zeros(np.shape(fftmap))
                std_phasemap_mean = np.zeros(np.shape(fftmap))
                # INV. ABEL TRANSF. MAP
                tgt_abelmap_mean = np.zeros(np.shape(fftmap))
                std_abelmap_mean = np.zeros(np.shape(fftmap))
                # PLASMA DENSITY
                tgt_dens_mean = np.zeros(np.shape(fftmap))
                std_dens_mean = np.zeros(np.shape(fftmap))

        '''
        BUILDING 2D AND 1D PLOTS
        '''
        # Plots are building from user select
        if values['fftradio'] == True:  # Plot FFT map result
            matrix_plot = fftmap
            indmax = np.zeros(2)
        elif values['filterradio'] == True:  # Plot gaussian filter map
            matrix_plot = gfilter
            indmax = np.zeros(2)
        elif values['phaseradio'] == True:  # Plot phase map result
            try:
                maxval = np.max(tgt_phasemap_mean)
                indmax = np.unravel_index(np.argmax(tgt_phasemap_mean, axis=None), tgt_phasemap_mean.shape)
                maxstd = std_phasemap_mean[indmax]
            except:
                maxval = 0
                maxstd = 0

            strmax = r'$\Delta\phi_{max}=%.3f \pm %.3f$' % (maxval, maxstd)

            if values['-checkstd-'] == False:
                matrix_plot = tgt_phasemap_mean
            else:
                matrix_plot = std_phasemap_mean
        elif values['abelradio'] == True:  # Plot gas density profile from IAT
            # max abs value and max abs position of phase
            try:
                maxval = np.min(tgt_abelmap_mean)
                indmax = np.unravel_index(np.argmax(tgt_abelmap_mean, axis=None), tgt_abelmap_mean.shape)
                maxstd = std_abelmap_mean[indmax]
            except:
                maxval = 0
                maxstd = 0

            strmax = r'$\Delta\phi_{r. max}=%.5f \pm %.5f$' % (maxval, maxstd)
            if values['-checkstd-'] == False:
                matrix_plot = tgt_abelmap_mean
            else:
                matrix_plot = std_abelmap_mean
        elif values['densradio'] == True:  # Plot gas density profile
            # max abs value and max abs position of phase
            try:
                maxval = np.max(tgt_dens_mean)
                indmax = np.unravel_index(np.argmax(tgt_dens_mean, axis=None), tgt_dens_mean.shape)
                maxstd = 100 * std_dens_mean[indmax] / maxval
            except:
                maxval = 0
                maxstd = 0
            strmax = r'$N_{g.max}=%.2e (\pm %.1f\%%)$' % (maxval, maxstd)

            if values['-checkstd-'] == False:
                matrix_plot = tgt_dens_mean
            else:
                matrix_plot = std_dens_mean

        rangeh, rangev = np.shape(matrix_plot)
        # Ajust slider for horizontal/vertical
        if values['-comboaxisymm-'] == 'vertical':
            window['sliderh'].update(range=(0, rangeh - 1), value=int(rangeh - indmax[0]))
        elif values['-comboaxisymm-'] == 'horizontal':
            window['sliderh'].update(range=(0, rangev - 1), value=indmax[1])

        newcmp = cm.get_cmap(values['-cmapcombo-'], 256)

        # Creating plot figure parameters
        fig, ax1 = plt.subplots(figsize=(4.9, 4))

        if values['filterradio'] == True:
            ax1.set_xlabel('$\\nu_x\hspace{.5}(arb.u.)$', fontsize=12)
            ax1.set_ylabel('$\\nu_y\hspace{.5}(arb.u.)$', fontsize=12)
            ax1.imshow(matrix_plot, cmap='gray')

        elif values['fftradio'] == True:
            ax1.set_xlabel('$\\nu_x\hspace{.5}(arb.u.)$', fontsize=12)
            ax1.set_ylabel('$\\nu_y\hspace{.5}(arb.u.)$', fontsize=12)
            ax1.imshow(matrix_plot, cmap='gray')
            if values['-combofringe-'] == 'horizontal':
                ax1.axhline(y=centerfilter, lw=1, alpha=0.5, color='red')
            else:
                ax1.axvline(x=centerfilter, lw=1, alpha=0.5, color='red')


        elif values['phaseradio'] == True:
            divider = make_axes_locatable(ax1)
            extentplot = np.shape(matrix_plot)
            x_max = extentplot[1] * factor * 1e6
            y_max = extentplot[0] * factor * 1e6
            cax = divider.append_axes("right", size="5%", pad=0.05)
            map_plot = ax1.imshow(matrix_plot, extent=[0, x_max, 0, y_max], cmap=newcmp)
            cb1 = fig.colorbar(map_plot, cax=cax)
            ax1.set_xlabel('$x\hspace{.5}(\mu m)$', fontsize=12)
            ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
            ax1.text(0, 1.05 * y_max, strmax)
            if values['-checkstd-'] == False:
                cb1.set_label(label='$\Delta\phi\hspace{.5} (rad)$', size=12, weight='bold')
            else:
                cb1.set_label(label='$\sigma_{\Delta\phi}\hspace{.5} (rad)$', size=12, weight='bold')
            if values['-comboaxisymm-'] == 'horizontal':
                ax1.axhline(y=axis_pos * factor * 1e6, lw=1, ls='--', alpha=0.5, color='black')
            else:
                ax1.axvline(x=axis_pos * factor * 1e6, lw=1, ls='--', alpha=0.5, color='black')

        elif values['densradio'] == True:
            divider = make_axes_locatable(ax1)
            extentplot = np.shape(matrix_plot)
            x_max = extentplot[1] * factor * 1e6
            y_max = extentplot[0] * factor * 1e6
            cax = divider.append_axes("right", size="5%", pad=0.05)
            map_plot = ax1.imshow(matrix_plot, extent=[-x_max/2, x_max/2, 0, y_max], cmap=newcmp)
            cb1 = fig.colorbar(map_plot, cax=cax)
            ax1.set_xlabel('$r\hspace{.5}(\mu m)$', fontsize=12)
            ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
            ax1.text(-x_max/2, 1.05 * y_max, strmax)
            if values['-checkstd-'] == False:
                cb1.set_label(label='$N_g\hspace{.5} (cm^{-3})$', size=12, weight='bold')
            else:
                cb1.set_label(label='$\sigma_{N_g}\hspace{.5} (cm^{-3})$', size=12, weight='bold')

        else:
            divider = make_axes_locatable(ax1)
            extentplot = np.shape(matrix_plot)
            x_max = extentplot[1] * factor * 1e6
            y_max = extentplot[0] * factor * 1e6
            cax = divider.append_axes("right", size="5%", pad=0.05)

            map_plot = ax1.imshow(matrix_plot, extent=[-x_max/2, x_max/2, 0, y_max], cmap=newcmp)
            cb1 = fig.colorbar(map_plot, cax=cax)
            ax1.set_xlabel('$r\hspace{.5}(\mu m)$', fontsize=12)
            ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
            ax1.text(-x_max/2, 1.05 * y_max, strmax)
            if values['-checkstd-'] == False:
                cb1.set_label(label='$\Delta\phi_{r}\hspace{.5} (rad/ \mu m)$', size=12, weight='bold')
            else:
                cb1.set_label(label='$\sigma_{\Delta\phi_{r}}\hspace{.5} (rad/ \mu m)$', size=12, weight='bold')

        fig.tight_layout(pad=2)
        fig_canvas_agg = draw_figure(window['canvasabel'].TKCanvas, fig)

        visible_f1d = False
        # Enable/Disable specific buttons and frames for 2D analysis
        window['Save Data'].update(disabled=False)
        window['Save Plot'].update(disabled=False)
        window['frame1d'].update(visible=False)
        window['2D Profile'].update(disabled=False)
        window['1D Profile'].update(disabled=False)
    #########################################################################
    # BUTTON DENS.PROFILE 2
    if event == '2D Profile':
        # Cleaning plots
        try:
            fig_canvas_agg.get_tk_widget().forget()
        except:
            plt.close('all')
        # set height position
        h_prof = -1.0

        # Plots are building from user select
        if values['fftradio'] == True:  # Plot FFT map result
            matrix_plot = fftmap
            indmax = np.zeros(2)
        elif values['filterradio'] == True:  # Plot gaussian filter map
            matrix_plot = gfilter
            indmax = np.zeros(2)
        elif values['phaseradio'] == True:  # Plot phase map result
            # max abs value and max abs position of phase
            try:
                maxval = np.max(tgt_phasemap_mean)
                indmax = np.unravel_index(np.argmax(tgt_phasemap_mean, axis=None), tgt_phasemap_mean.shape)
                maxstd = std_phasemap_mean[indmax]
            except:
                maxval = 0
                maxstd = 0
            strmax = r'$\Delta\phi_{max}=%.3f \pm %.3f$' % (maxval, maxstd)

            if values['-checkstd-'] == False:
                matrix_plot = tgt_phasemap_mean
            else:
                matrix_plot = std_phasemap_mean
        elif values['abelradio'] == True:  # Plot gas density profile from IAT
            # max abs value and max abs position of phase
            try:
                maxval = np.max(tgt_abelmap_mean)
                indmax = np.unravel_index(np.argmax(tgt_abelmap_mean, axis=None), tgt_abelmap_mean.shape)
                maxstd = std_abelmap_mean[indmax]
            except:
                maxval = 0
                maxstd = 0
            strmax = r'$\Delta\phi_{r. max}=%.5f \pm %.5f$' % (maxval, maxstd)
            if values['-checkstd-'] == False:
                matrix_plot = tgt_abelmap_mean
            else:
                matrix_plot = std_abelmap_mean
        elif values['densradio'] == True:  # Plot gas density profile
            # max value and max value position of density
            try:
                maxval = np.max(tgt_dens_mean)
                indmax = np.unravel_index(np.argmax(tgt_dens_mean, axis=None), tgt_dens_mean.shape)
                maxstd = 100 * std_dens_mean[indmax] / maxval
            except:
                maxval = 0
                maxstd = 0
            strmax = r'$N_{g.max}=%.2e (\pm %.1f\%%)$' % (maxval, maxstd)

            if values['-checkstd-'] == False:
                matrix_plot = tgt_dens_mean
            else:
                matrix_plot = std_dens_mean

        # colormap distribution
        newcmp = cm.get_cmap(values['-cmapcombo-'], 256)

        try:
            # clearing figures and plots
            fig_canvas_agg.get_tk_widget().forget()
            plt.close('all')

            # Instead of plt.show
            fig, ax1 = plt.subplots(figsize=(4.9, 4))

            if values['filterradio'] == True:
                ax1.set_xlabel('$\\nu_x\hspace{.5}(arb.u.)$', fontsize=12)
                ax1.set_ylabel('$\\nu_y\hspace{.5}(arb.u.)$', fontsize=12)
                abel_plot = ax1.imshow(matrix_plot, cmap='gray')

            elif values['fftradio'] == True:
                ax1.set_xlabel('$\\nu_x\hspace{.5}(arb.u.)$', fontsize=12)
                ax1.set_ylabel('$\\nu_y\hspace{.5}(arb.u.)$', fontsize=12)
                if values['-combofringe-'] == 'horizontal':
                    ax1.axhline(y=centerfilter, lw=1, alpha=0.5, color='red')
                else:
                    ax1.axvline(x=centerfilter, lw=1, alpha=0.5, color='red')
                ax1.imshow(matrix_plot, cmap='gray')

            elif values['phaseradio'] == True:
                divider = make_axes_locatable(ax1)
                extentplot = np.shape(matrix_plot)
                x_max = extentplot[1] * factor * 1e6
                y_max = extentplot[0] * factor * 1e6
                cax = divider.append_axes("right", size="5%", pad=0.05)
                abel_plot = ax1.imshow(matrix_plot, extent=[0, x_max, 0, y_max], cmap=newcmp)
                cb1 = fig.colorbar(abel_plot, cax=cax)
                ax1.set_xlabel('$x\hspace{.5}(\mu m)$', fontsize=12)
                ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
                ax1.text(0, 1.1 * y_max, strmax)
                if values['-checkstd-'] == False:
                    cb1.set_label(label='$\Delta\phi\hspace{.5} (rad)$', size=12, weight='bold')
                else:
                    cb1.set_label(label='$\sigma_{\Delta\phi}\hspace{.5} (rad)$', size=12, weight='bold')

            elif values['densradio'] == True:
                divider = make_axes_locatable(ax1)
                extentplot = np.shape(matrix_plot)
                x_max = extentplot[1] * factor * 1e6
                y_max = extentplot[0] * factor * 1e6
                cax = divider.append_axes("right", size="5%", pad=0.05)
                abel_plot = ax1.imshow(matrix_plot, extent=[-x_max/2, x_max/2, 0, y_max], cmap=newcmp)
                cb1 = fig.colorbar(abel_plot, cax=cax)
                ax1.set_xlabel('$r\hspace{.5}(\mu m)$', fontsize=12)
                ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
                ax1.text(-x_max/2, 1.1 * y_max, strmax)
                if values['-checkstd-'] == False:
                    cb1.set_label(label='$N\hspace{.5} (cm^{-3})$', size=12, weight='bold')
                else:
                    cb1.set_label(label='$\sigma_{N}\hspace{.5} (cm^{-3})$', size=12, weight='bold')

            else:
                divider = make_axes_locatable(ax1)
                extentplot = np.shape(matrix_plot)
                x_max = extentplot[1] * factor * 1e6
                y_max = extentplot[0] * factor * 1e6
                cax = divider.append_axes("right", size="5%", pad=0.05)

                abel_plot = ax1.imshow(matrix_plot, extent=[-x_max/2, x_max/2, 0, y_max], cmap=newcmp)
                cb1 = fig.colorbar(abel_plot, cax=cax)
                ax1.set_xlabel('$r\hspace{.5}(\mu m)$', fontsize=12)
                ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
                ax1.text(-x_max/2, 1.1 * y_max, strmax)
                if values['-checkstd-'] == False:
                    cb1.set_label(label='$\Delta\phi_{r}\hspace{.5} (rad/ \mu m)$', size=12, weight='bold')
                else:
                    cb1.set_label(label='$\sigma_{\Delta\phi_{r}}\hspace{.5} (rad/ \mu m)$', size=12, weight='bold')

            fig.tight_layout(pad=2)
            fig_canvas_agg = draw_figure(window['canvasabel'].TKCanvas, fig)

            visible_f1d = False
            window['frame1d'].update(visible=False)

        except:
            continue
    #########################################################################
    # BUTTON DENS.PROFILE 1D AND SLIDER POSITION
    if (event == '1D Profile') or (event == 'sliderh'):
        # Cleaning plots
        try:
            fig_canvas_agg.get_tk_widget().forget()
        except:
            plt.close('all')
        # set height position
        h_prof = -1.0

        if values['-checkstd-'] == True:
            window['-checkpos1-'].update(False)
            window['-checkpos2-'].update(False)
            window['-checkpos3-'].update(False)
        # Plots are building from user select
        if values['fftradio'] == True:  # Plot FFT map result
            matrix_plot = fftmap
            matrix_plot_std = np.zeros(np.shape(matrix_plot))
            headerfile = 'pixel, R.Int.(arb.u)'
        elif values['filterradio'] == True:  # Plot gaussian filter map
            matrix_plot = gfilter
            matrix_plot_std = np.zeros(np.shape(matrix_plot))
            headerfile = 'pixel, R.Int.(arb.u)'
        elif values['phaseradio'] == True:  # Plot phase map result
            matrix_plot = tgt_phasemap_mean
            matrix_plot_std = std_phasemap_mean
            headerfile = '\nPositions(µm), acc. phase(rad)'
        elif values['abelradio'] == True:  # Plot gas density profile from IAT
            matrix_plot = tgt_abelmap_mean
            matrix_plot_std = norm_phasemap
            headerfile = '\nPositions(µm), radial phase(rad/m)'
        elif values['densradio'] == True:  # Plot gas density profile
            headerfile = '\nPosition(µm), Ng (1/cm³)'
            matrix_plot = tgt_dens_mean
            matrix_plot_std = std_dens_mean

        rangeh, rangev = np.shape(matrix_plot)

        try:
            # clearing figures and plots
            fig_canvas_agg.get_tk_widget().forget()
            plt.close('all')

            # create r axis according to symmetry in micrometers (µm)
            if values['-comboaxisymm-'] == 'vertical' or 'none':
                raxis = np.arange(-rangev / 2, rangev / 2, 1)
                raxis_um = raxis * factor * 1e6  # um
                # set origin position (exit nozzle position) and slider position
                h_prof = 0
                pos_0 = rangeh - 1
                pos = pos_0 - int(values['sliderh'])
                array_plot = matrix_plot[pos]
                array_std = matrix_plot_std[pos]

            elif values['-comboaxisymm-'] == 'horizontal':
                raxis = np.arange(-rangeh / 2, rangeh / 2, 1)
                raxis_um = raxis * factor * 1e6  # um
                # set origin position (exit nozzle position) and slider position
                h_prof = 0
                pos_0 = 0
                pos = int(values['sliderh'])
                array_plot = matrix_plot[:, pos]
                array_std = matrix_plot_std[:, pos]

            # convert vertical array positions to height positions in µm
            h_prof = int(values['sliderh']) * factor * 1e6

            h_prof1 = int(get_value('-pos1-', values))
            pos1 = int(h_prof1 / (factor * 1e6))
            h_prof2 = int(get_value('-pos2-', values))
            pos2 = int(h_prof2 / (factor * 1e6))
            h_prof3 = int(get_value('-pos3-', values))
            pos3 = int(h_prof3 / (factor * 1e6))

            # Creating plot parameters
            fig, ax1 = plt.subplots(figsize=(4.9, 4))

            ax1.set_xlabel('$r\hspace{.5}(\mu m)$', fontsize=12)

            if values['densradio'] == True:
                labelplot = '$%d \hspace{.5}\mu m$'
                ax1.plot(raxis_um, array_plot, label=labelplot % h_prof, lw=2, color="blue")
                ax1.set_ylabel('$N_g\hspace{.5} (cm^{-3})$', fontsize=12)
                strmax = r'$N_{g.max}=%.2e (\pm %.1f\%%)$' % (np.max(array_plot), 100 * array_std[np.argmax(array_plot)] / np.max(array_plot))
                ax1.text(raxis_um[0] / 2, 1.1 * np.max(array_plot), strmax)
                ax1.set_ylim(bottom=0.)
                try:
                    if values['-checkpos1-'] == False and values['-checkpos2-'] == False and values['-checkpos3-'] == False:
                        FWHM, yFWHM, x1, x2 = peak_widths(np.abs(array_plot), [np.argmax(array_plot)], rel_height=0.5)
                        ax1.hlines(yFWHM, raxis_um[int(x1)], raxis_um[int(x2)], label='%.1f $\mu m (FWHM)$'%(FWHM[0]*factor*1e6), lw=1, color="black", ls='--')
                except:
                    continue
                if values['-checkstd-'] == True:
                    ax1.errorbar(raxis_um, array_plot, yerr=array_std, label='$\sigma_{N_g}$', alpha=0.2,
                                 color="blue")

            if values['abelradio'] == True:
                labelplot = '$(\Delta\phi_r)_{%d \hspace{.5}\mu m}$'
                ax1.plot(raxis_um, array_plot, label=labelplot % h_prof, lw=2, color="blue")
                ax1.set_ylabel('$\Delta\phi_{r}\hspace{.5} (rad/ \mu m)$', fontsize=12)
                ax1.set_ylim(bottom=0.)
                if values['-checkstd-'] == True:
                    ax1.plot(raxis_um, array_std, '--', label='$(\Delta\phi_{norm})_{%d \hspace{.5}\mu m}$' % h_prof,
                             lw=2, color="red")
                    ax1.fill_between(raxis_um, array_plot, array_std, color="orange", alpha=0.5,
                                     label='$(\sigma_{Abel})_{%d}$' % h_prof)

            if values['phaseradio'] == True:
                labelplot = '$%d \hspace{.5}\mu m$'
                ax1.plot(raxis_um, array_plot, label=labelplot % h_prof, lw=2, color="blue")
                ax1.set_ylabel('$\Delta\phi\hspace{.5} (rad)$', fontsize=12)
                ax1.set_ylim(bottom=0.)
                if values['-checkstd-'] == True:
                    ax1.errorbar(raxis_um, array_plot, yerr=array_std, label='$\sigma_{\Delta\phi}$', alpha=0.2,
                                 color="blue")

            if values['filterradio'] == True or values['fftradio'] == True:
                if values['-combofringe-'] == 'vertical':
                    ax1.plot(np.arange(0, rangev), summap, label='FFT Profile', lw=2, color="gray")
                    if values['filterradio'] == True:
                        ax1.axvline(x=centerfilter + sigma_gfilter, label='$\Delta\\nu$', lw=1, alpha=0.5, color='red')
                        ax1.axvline(x=centerfilter - sigma_gfilter, lw=1, alpha=0.5, color='red')
                    else:
                        ax1.scatter(filterpoints, [summap[i] for i in filterpoints], color='red', marker='^',
                                    label='$\\nu_0$')
                        [ax1.text(x, summap[x], '%d' % x) for x in filterpoints]

                else:
                    ax1.plot(np.arange(0, rangeh), summap, label='FFT Profile', color="gray")
                    if values['filterradio'] == True:
                        ax1.axvline(x=centerfilter + sigma_gfilter, label='$\Delta\\nu$', lw=1, alpha=0.5,
                                    color='red')
                        ax1.axvline(x=centerfilter - sigma_gfilter, lw=1, alpha=0.5, color='red')
                    else:
                        ax1.scatter(filterpoints, [summap[i] for i in filterpoints], color='red', marker='^',
                                    label='$\\nu_0$')
                        [ax1.text(x, summap[x], '%d' % x) for x in filterpoints]
                ax1.set_xlabel('$Frequency\hspace{.5}\\nu\hspace{.5}(pixel)$', fontsize=12)
                ax1.set_ylabel('$Relative\hspace{.5}Intensity\hspace{.5} (a.u.)$', fontsize=12)


            # Including new 1D density profile for another height from origin height position
            if values['-checkpos1-'] == True and values['-checkstd-'] == False:
                h_prof1 = int(get_value('-pos1-', values))
                pos1 = int(h_prof1 / (factor * 1e6))
                if values['-comboaxisymm-'] == 'vertical':
                    ax1.plot(raxis_um, matrix_plot[pos_0 - pos1], label=labelplot % (h_prof1), lw=1,
                             color="purple")
                elif values['-comboaxisymm-'] == 'horizontal':
                    ax1.plot(raxis_um, matrix_plot[:, pos1], label=labelplot % (h_prof1), lw=1,
                             color="purple")

            if values['-checkpos2-'] == True and values['-checkstd-'] == False:
                h_prof2 = int(get_value('-pos2-', values))
                pos2 = int(h_prof2 / (factor * 1e6))
                if values['-comboaxisymm-'] == 'vertical':
                    ax1.plot(raxis_um, matrix_plot[pos_0 - pos2], label=labelplot % (h_prof2), lw=1,
                             color="red")
                elif values['-comboaxisymm-'] == 'horizontal':
                    ax1.plot(raxis_um, matrix_plot[:, pos2], label=labelplot % h_prof2, lw=1,
                             color="red")

            if values['-checkpos3-'] == True and values['-checkstd-'] == False:
                h_prof3 = int(get_value('-pos3-', values))
                pos3 = int(h_prof3 / (factor * 1e6))
                if values['-comboaxisymm-'] == 'vertical':
                    ax1.plot(raxis_um, matrix_plot[pos_0 - pos3], label=labelplot % (h_prof3), lw=1,
                             color="orange")
                elif values['-comboaxisymm-'] == 'horizontal':
                    ax1.plot(raxis_um, matrix_plot[:, pos3], label=labelplot % h_prof3, lw=1,
                             color="orange")

            ax1.legend()
            ax1.grid(True)
            fig.tight_layout(pad=2)
            fig_canvas_agg = draw_figure(window['canvasabel'].TKCanvas, fig)
            visible_f1d = True
            window['frame1d'].update(visible=visible_f1d)
        except:
            continue
    #########################################################################
    # Saving results
    #########################################################################
    #  BUTTON SAVEPLOT
    elif event == 'Save Plot':
        save_filename_plot = sg.popup_get_file('File',
                                               file_types=[("PNG (*.png)", "*.png"), ("All files (*.*)", "*.*")],
                                               save_as=True, no_window=True)
        if save_filename_plot:
            # save the plot
            plt.savefig(save_filename_plot)
            sg.popup(f"Saved: {save_filename_plot}")

    ########################################################################
    #  BUTTON SAVE DATA
    elif event == 'Save Data':
        save_filename_data = sg.popup_get_file('File',
                                               file_types=[("DAT (*.dat)", "*.dat"), ("TXT (*.txt)", "*.txt")],
                                               save_as=True, no_window=True)

        if save_filename_data:
            file_data = open(save_filename_data, 'a')
            file_data.seek(0)  # sets  point at the beginning of the file
            file_data.truncate()

            # Saving 1D plots
            if visible_f1d == True:
                # save data plot
                file_data.write(headerfile)
                if (values['filterradio'] == True or values['fftradio'] == True):
                    if values['-combofringe-'] == 'vertical':
                        list_data = np.vstack((np.arange(0, nrmap, 1), summap))
                    elif values['-combofringe-'] == 'horizontal':
                        # sum of array lines (horizontal)
                        list_data = np.vstack((np.arange(0, nlmap, 1), summap))

                if h_prof >= 0.0 and (
                        values['phaseradio'] == True or values['abelradio'] == True or values['densradio'] == True):
                    file_data.write(',for height(s) on axisymmetric of %.0f (µm)' % h_prof)
                    list_data = np.vstack((raxis_um, array_plot))
                    if values['-checkstd-'] == True:
                        file_data.write('with standard deviation')
                        list_data = np.vstack((list_data, array_std))
                    # verify additional height positions 1, 2 and 3
                    else:
                        if values['-checkpos1-'] == True:
                            if values['-comboaxisymm-'] == 'vertical':
                                list_data = np.vstack((list_data, matrix_plot[pos_0 - pos1]))
                            else:
                                list_data = np.vstack((list_data, matrix_plot[:, pos1]))
                            file_data.write(',%.0f (µm)' % h_prof1)

                        if values['-checkpos2-'] == True:
                            if values['-comboaxisymm-'] == 'vertical':
                                list_data = np.vstack((list_data, matrix_plot[pos_0 - pos2]))
                            else:
                                list_data = np.vstack((list_data, matrix_plot[:, pos2]))
                            file_data.write(',%.0f (µm)' % h_prof2)

                        if values['-checkpos3-'] == True:
                            if values['-comboaxisymm-'] == 'vertical':
                                list_data = np.vstack((list_data, matrix_plot[pos_0 - pos3]))
                            else:
                                list_data = np.vstack((list_data, matrix_plot[:, pos3]))
                            file_data.write(',%.0f (µm)' % h_prof3)

                file_data.write('\n')
                list_str = (np.transpose(list_data))
                np.savetxt(file_data, list_str, fmt='%.2e', delimiter=',')
                file_data.close()

            # Save 2D array
            else:
                if values['-checkstd-'] == True:
                    np.savetxt(file_data, matrix_plot_std, fmt='%.3e', delimiter=',')
                else:
                    np.savetxt(file_data, matrix_plot, fmt='%.3e', delimiter=',')

                sg.popup(f"Saved: {save_filename_data}")
                file_data.close()
        else:
            continue

window.close()

