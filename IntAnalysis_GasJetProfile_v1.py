# Software: Interfetometry Analysis - Gas-Jet Profile (Version 1.0)
# Authors: Jhonatha Ricardo dos Santos, Armando Zuffi, Ricardo Edgul Samad, Nilson Dias Vieira Junior
# Python 3.11

# LYBRARIES
# Numpy from https://numpy.org/
# PyAbel/PyAbel:v0.9.0rc1 from https://doi.org/10.5281/zenodo.7401589.svg
# PySimpleGUI from pysimplegui.org
# Matplotlib from matplotlib.org
# Scipy from scipy.org
# Scikit-image from  https://doi.org/10.7717/peerj.453
# Pillow (PIL Fork) 9.3.0 from https://pypi.org/project/Pillow
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

from io import BytesIO
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import cm
from matplotlib.colors import ListedColormap
from scipy.ndimage import gaussian_filter
from scipy.signal import peak_widths, find_peaks
from skimage.restoration import unwrap_phase
from PIL import Image, ImageDraw, UnidentifiedImageError

# Matplotlib Tk style
matplotlib.use('TkAgg')
# Font and theme of PysimpleGUI
AppFont = 'Arial 16 bold'
sg.theme('DarkGrey4')

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
polargas = '1.710'  #
sigma_gfilter = '10'  # sigma of gauss function
centerfilter = '0'  # Center of the gaussian filter application
specificheat = '1.47'  # specific heat of gas in A³ (Angstrom)
factor = '1.000'  # factor um/pixel
# Image parameters
h_prof = -1.0  # heigth null
rotate_degree = 0.0  # angle to image rotation
# Initial values to cut image
begin_x = '100'
begin_y = '100'
end_x = '300'
end_y = '300'
# Initial values of heigths for 1D analysis
pos1 = '50'
pos2 = '75'
pos3 = '100'

# Images Dimensions
width, height = size = 428, 342  # Scale image - interferogram
width2, height2 = size2 = 214, 172  # Scale image - Ref
width3, height3 = size3 = 428, 342  # Scale image - Result
# Min and Max values of Interferogram Image
minvalue_x, maxvalue_x, minvalue_y, maxvalue_y = 0, 428, 0, 342
#Frame 1D visible
visible_f1d=False

#################################################################################
# FUNCTIONS
#################################################################################
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
# DRAW SQUARE IN CANVAS
def draw_figure(canvas, figure):
    '''
    Drawing rectangle figure on canvas
    :param canvas:
    :param figure:
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

    if os.path.exists(image_file):
        shutil.copy(image_file, tmp_file)
        imagetmp = Image.open(tmp_file)
        imagetmp = imagetmp.resize(size)
        draw = ImageDraw.Draw(imagetmp)
        draw.rectangle((begin_x, begin_y, end_x, end_y), width=2,
                       outline='#FFFFFF')  ##DCDCDC
        imagetmp.save(tmp_file)
        bio = io.BytesIO()
        imagetmp.save(bio, format='PNG')
        window["image1"].update(data=bio.getvalue(), size=size)
    # NEW COLORMAPS
# COLORMAP DEFINITION
def func_colormap(n):
    '''
    Color distribution at the colormap
    :param n: n_order for colormap (linear,quadratic, cubic)
    :return:
    '''
    return ListedColormap(cm.get_cmap('rainbow_r', 512)(np.power(np.linspace(1, 0, 512), n)))
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
# CREATING INTENSITY DIST OF THE FRINGES
def intensity_dist(data, fringe_axis):
    '''
    Calculate the intensity distribution of 2D array
    :param: 2D arrays and arientation of fringes
    :return: 2D array
    '''
    data_dist = np.zeros(np.shape(data))
    nl, nr = np.shape(data)
    # vertical fringes
    if fringe_axis == 0:
        for i in range(0, nl):
            y = data[i, :]
            x = np.arange(len(y))
            ypeaks, _ = find_peaks(y)
            data_dist[i, :] = np.interp(x, ypeaks, y[ypeaks])
            # horizontal fringes
    elif fringe_axis == 1:
        for i in range(0, nr):
            y = data[:, i]
            x = np.arange(len(y))
            ypeaks, _ = find_peaks(y)
            data_dist[:, i] = np.interp(x, ypeaks, y[ypeaks])
    return data_dist
# CREATING SHIFT AND WIDTHS OF THE FRINGES
def fringes_info(data1, data2, data3):
    '''
     Calculate 2D array shifts and widths fringes distribution
    :param n: slice 2D arrays, slice 2D array of ref. image, 2D array ref.
    :return: 2D arrays of shifts and widths fringes distribution
    '''
    data_shift = np.zeros(np.shape(data3))
    data_dist = np.zeros(np.shape(data3))
    nl, nr = np.shape(data3)
    ypeaks2, _ = find_peaks(data2)
    ypeaks1, _ = find_peaks(data1)
    for i in range(0, len(ypeaks1)):
        try:
            teste = np.isclose(ypeaks2[:i + 1], ypeaks1[:i + 1], rtol=0, atol=2.1)
            if teste[i] == False and ypeaks1[i] > ypeaks2[i]:
                ypeaks2 = np.delete(ypeaks2, i)
                i = 0
            if teste[i] == False and ypeaks1[i] < ypeaks2[i]:
                ypeaks1 = np.delete(ypeaks1, i)
                i = 0
        except:
            break
    while len(ypeaks1) > len(ypeaks2):
        ypeaks1 = np.delete(ypeaks1, -1)
    while len(ypeaks1) < len(ypeaks2):
        ypeaks2 = np.delete(ypeaks2, -1)
    x = np.arange(nr)
    dist_i = np.interp(x, np.arange(len(np.diff(ypeaks1))), np.diff(ypeaks1))
    shift_i = np.interp(x, np.arange(len(ypeaks1)), (abs(ypeaks1 - ypeaks2) - np.min(abs(ypeaks1 - ypeaks2))))
    for j in range(0, nl):
        data_shift[j] = shift_i
        data_dist[j] = dist_i
    return data_shift, data_dist

'''
########################################################################################
#WINDOWS LAYOUT
Building frames for main windows
########################################################################################
'''
# LAYOUT INTERFEROMETER IMAGE
layout_frame_ImgSample = [
    [sg.Image(size=size, background_color='black',
              key='image1', expand_x=True)
     ],
    [sg.Input(expand_x=True, disabled=True, key='file1', visible='True')],
    [sg.Button('Open File', font='Arial 10 bold'),
     sg.Button('Rotate (°)', visible=True, font='Arial 10 bold', disabled=True),
     sg.Input('0', size=(5, 1), key='-DEGREE-', enable_events=True),
     sg.Text('  Image Rescale (w,h):'),
     sg.Text('', key='-scale1-')],
]
# LAYOUT REFERENCE IMAGE
layout_frame_ImgReference = [
    [sg.Image(size=size2, background_color='black',
              key='image2', enable_events=True)],
    [sg.Input(expand_x=True, disabled=True, key='file2', visible='True')],
    [sg.Button('Open Ref.', font='Arial 10 bold')],
]
# LAYOUT SELECT COORD. AREA OPTIONS
layout_area_coord = [
    [sg.Text('Coord X'),
     sg.Spin([i for i in range(minvalue_x, maxvalue_x + 1)], initial_value=begin_x, key='-BEGIN_X-', size=(4, 1),
             enable_events=True),
     sg.Spin([i for i in range(minvalue_x, maxvalue_x + 1)], initial_value=end_x, key='-END_X-', size=(4, 1),
             enable_events=True)],

    [sg.Text('Coord Y'),
     sg.Spin([i for i in range(minvalue_y, maxvalue_y + 1)], initial_value=begin_y, key='-BEGIN_Y-', size=(4, 1),
             enable_events=True),
     sg.Spin([i for i in range(minvalue_y, maxvalue_y + 1)], initial_value=end_y, key='-END_Y-', size=(4, 1),
             enable_events=True)],

]
# LAYOUT SELECT AREA OPTIONS
layout_area_selection = [
    [sg.Button('Select Analysis Area', size=(20, 2), font='Arial 10 bold', disabled=True)],
    [sg.Checkbox('Cut selected area', default=True, key='-checkcut-')],
    [sg.Frame('Area Coord.', layout_area_coord, size=(178, 80), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 10 bold')],
]
# LAYOUT INPUT GAS AND RADIATION PARAMETERS
layout_input_parameters = [
    [sg.Text('Laser \nWavelength (nm): '),
     sg.Input(lambda0, size=(5, 1), key='-lambda0-', enable_events=True)],
    [sg.Text('Unc. laser\nWavelength (nm): '),
     sg.Input(unc_lambda0, size=(5, 1), key='-unclambda0-', enable_events=True)],
    [sg.Text('Gas Type:           '),
     sg.Combo(['H2', 'N2', 'He', 'Ar', '--'], default_value='N2', key='-combogas-', enable_events=True)],
    [sg.Text('Polarizability (Å³):'),
     sg.Input(polargas, size=(5, 1), key='-polargas-', enable_events=True)],
    [sg.Text('Specific Heat:     '),
     sg.Input(specificheat, size=(5, 1), key='-specificheatgas-', enable_events=True)],
]
# LAYOUT INPUT MEASUREMENT PARAMETERS
layout_analysis_parameters = [
    [sg.Text('Scaling Factor (µm/pixel):           '),
     sg.Input(factor, size=(5, 1), key='-factor-', enable_events=True)],
    [sg.Text('Sigma - Gaussian Blur (pixel):     '),
     sg.Input(sigma_gfilter, size=(5, 1), key='-sigma_gfilter-', enable_events=True)],
    [sg.Text('Gaussian Filter position (pixel):   '),
     sg.Input(centerfilter, size=(5, 1), key='-centerfilter-', enable_events=True)],
    [sg.Text('Fringes Orientation:       '),
     sg.Combo(['vertical', 'horizontal'], default_value='vertical', key='-combofringe-')],
    [sg.Text('Axisymmetric:              '),
     sg.Combo(['vertical', 'horizontal'], default_value='vertical', key='-comboaxisymm-')]
]
# LAYOUT FRAME OF ALL INPUT OPTIONS
layout_frame_Options = [
    [sg.Frame('Select Area', layout_area_selection, size=(198, 210), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 10 bold'),
     sg.Frame('Input Parameters', layout_input_parameters, size=(178, 210), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 10 bold'),
     sg.Frame('Analysis Parameters', layout_analysis_parameters, size=(268, 210), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 10 bold')],
]
# LAYOUT FRAME LEFT - INPUTS
layout_frame_ImagesL = [
    [sg.Frame("Interferogram", layout_frame_ImgSample, size=(440, 430), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 12 bold')],
]
# LAYOUT FRAME RIGHT - INPUTS AND APPLY
layout_frame_ImagesR = [
    [sg.Frame("Reference", layout_frame_ImgReference, size=(258, 258), title_location=sg.TITLE_LOCATION_TOP,
              vertical_alignment="top", font='Arial 12 bold')],
    [sg.Button('Apply Algorithm', size=(30, 6), font='Arial 12 bold', disabled=True, button_color='gray')],
    [sg.Button('Exit', size=(30, 2), button_color='black', font='Arial 10 bold')]
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
    [sg.Checkbox('Prof. 1 (um)', default=False, key='-checkpos1-'),
     sg.Input(pos1, size=(4, 1), key='-pos1-', enable_events=True),
     sg.Checkbox('Prof. 2 (um)', default=False, key='-checkpos2-'),
     sg.Input(pos2, size=(4, 1), key='-pos2-', enable_events=True),
     sg.Checkbox('Prof. 3 (um)', default=False, key='-checkpos3-'),
     sg.Input(pos3, size=(4, 1), key='-pos3-', enable_events=True)],
    [sg.Slider(key='sliderh', range=(490, 0), orientation='h', size=(400, 20), default_value=0,
               enable_events=True)]
]
# LAYOUT STAGES OF THE TREATMENT
layout_frame_Steps = [
    [sg.Radio('Fourier \nTransform', "RADIO1", default=False, key='fftradio'),
     sg.Radio('Gaussian \nFilter', 'RADIO1', default=False, key='filterradio'),
     sg.Radio('Accumulated \nPhase', "RADIO1", default=True, key='phaseradio'),
     sg.Radio('Abel \nTransform', "RADIO1", default=False, key='abelradio'),
     sg.Radio('Density \nProfile', "RADIO1", default=False, key='densradio')],
]
# LAYOUT GLOBAL OUTPUTS
layout_frame_Result = [
    [sg.Frame('Stages', layout_frame_Steps, size=(490, 70), title_location=sg.TITLE_LOCATION_TOP_LEFT,
              key='framesteps', font='Arial 10 bold')],
    [sg.Button('Dens. Profile 1D', size=(16, 1), font='Arial 10 bold', disabled=True),
     sg.Button('Dens. Profile 2D', size=(16, 1), font='Arial 10 bold', disabled=True),
     sg.Checkbox('Standard deviation', default=False, key='-checkstd-'),
     ],
    [sg.Canvas(key='canvasabel', size=(490, 400), background_color='black')],
    [sg.Button('Save Plot', size=(16, 2), disabled=True, font='Arial 10 bold'),
     sg.Button('Save Data', size=(16, 2), disabled=True, font='Arial 10 bold'),
     sg.Text('Colormap Dist.:'),
     sg.Combo(['Linear', 'Quadratic', 'Cubic'], default_value='Linear', key='-cmapcombo-', enable_events=True)
     ],
    [sg.Frame('Density Profile - 1D (Axisymmetry)', layout_frame_plot1D, size=(490, 100), title_location=sg.TITLE_LOCATION_TOP_LEFT,
              visible=False, key='frame1d', font='Arial 10')],

]
# lAYOUT GLOBAL
layout = [
    [sg.Frame("Interferogram Images", layout_frame_Images, size=(700, 700), font='Arial 12 bold'),
     sg.Frame("Density Profile", layout_frame_Result, size=(500, 700), title_location=sg.TITLE_LOCATION_TOP,
              font='Arial 12 bold')],
]
######################################################################################################################
window = sg.Window("Interferogram Analysis - Gas-Jet Profile (Version 1.0)", layout, margins=(1, 1), finalize=True)
######################################################################################################################
'''
####################################################################################################
#Windows events
####################################################################################################
'''
while True:
    event, values = window.read()
    # Removing temp files when the main window is closed
    if event =='Exit' or event == sg.WINDOW_CLOSED:
        if '_temp.png' in path1:
            os.remove(path1)
            os.remove(path2)
        break
    ########################################################################
    # OPEN INTERFEROGRAM IMAGE
    elif event == "Open File":
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
                originalgassnp = Image.new(mode='L', size=(720, 576))
                originalgassnp.putdata(data0)

                ipath = ipath.replace('.snp', '_temp.png')

                originalgaspng = originalgassnp.save(ipath)
                if len(path_files) == 0:
                    path1 = ipath
                    window['file1'].update(path1)

                path_files.append(ipath)

        try:
            # Open Files
            originalgas = []
            for i in range(0, len(path_files)):
                originalgas.append(Image.open(path_files[i]))
            apply_drawing(values, window)

        except UnidentifiedImageError:
            continue

        # scale 1: scale for interferogram image
        w, h = originalgas[0].size
        scale = (round(width / w, 2), round(height / h, 2))

        if scale != 1:
            im1 = originalgas[0].resize(size)
        else:
            im1 = originalgas[0]
        data1 = image_to_data(im1)

        window['image1'].update(data=data1, size=size)
        window['-scale1-'].update(scale)

        # Enable buttons
        if path1 != '' and path2 != '':
            window['Rotate (°)'].update(disabled=False)
            window['-DEGREE-'].update(visible=True)
            window['Apply Algorithm'].update(disabled=False)
            window['Select Analysis Area'].update(disabled=False)
    ########################################################################
    # OPEN REFERENCE FILE
    elif event == "Open Ref.":
        path_file2 = sg.popup_get_file("", no_window=True)
        if path_file2:
            path2 = path_file2
        else:
            path2 = ''
        window['file2'].update(path2)
        # No file
        if path2 == '':
            continue
        # Create PNG files from files SNP
        if '.snp' in path2:
            databinary2 = getBinaryData(path2)
            data02 = np.flip(databinary2[60:60 + 720 * 576])
            originalrefsnp = Image.new(mode='L', size=(720, 576))
            originalrefsnp.putdata(data02)

            path2 = path2.replace('.snp', '_temp.png')

            originalrefpng = originalrefsnp.save(path2)

            window['file2'].update(path2)

        # Open files
        try:
            originalref = Image.open(path2)

        except UnidentifiedImageError:
            continue
        w2, h2 = originalref.size
        scale2 = (round(width2 / w2, 2), round(height2 / h2, 2))

        if scale2 != 1:
            im2 = originalref.resize(size2)
        else:
            im2 = originalref
        data2 = image_to_data(im2)
        window['image2'].update(data=data2, size=size2)
        # Enable buttons
        if path1 != '' and path2 != '':
            window['Rotate (°)'].update(disabled=False)
            window['-DEGREE-'].update(visible=True)
            window['Apply Algorithm'].update(disabled=False)
            window['Select Analysis Area'].update(disabled=False)
    ########################################################################
    # COMBO GAS TYPE
    elif event == '-combogas-':
        if values['-combogas-'] == 'N2':
            window['-polargas-'].update(value='1.710')  # cm3
            window['-specificheatgas-'].update(value='1.47')

        elif values['-combogas-'] == 'H2':
            window['-polargas-'].update(value='0.787')  # cm3
            window['-specificheatgas-'].update(value='1.405')

        elif values['-combogas-'] == 'He':
            window['-polargas-'].update(value='0.208')  # cm3
            window['-specificheatgas-'].update(value='1.667')

        elif values['-combogas-'] == 'Ar':
            window['-polargas-'].update(value='1.664')  # cm3
            window['-specificheatgas-'].update(value='1.667')
    ########################################################################
    # BUTTON SELECT AREA
    elif event == 'Select Analysis Area':
        apply_drawing(values, window)
        centerfilter = 0
        window['-centerfilter-'].update('0')
    # BUTTON COORD AREA
    elif event == '-BEGIN_X-':
        apply_drawing(values, window)
        centerfilter = 0
        window['-centerfilter-'].update('0')
    elif event == '-END_X-':
        apply_drawing(values, window)
        centerfilter = 0
        window['-centerfilter-'].update('0')
    elif event == '-BEGIN_Y-':
        apply_drawing(values, window)
        centerfilter = 0
        window['-centerfilter-'].update('0')
    elif event == '-END_Y-':
        apply_drawing(values, window)
        centerfilter = 0
        window['-centerfilter-'].update('0')
    #########################################################################
    # BUTTON ROTATE
    elif event == 'Rotate (°)':
        image_file = values["file1"]
        rotate_degree = (get_value('-DEGREE-', values))
        if os.path.exists(image_file):
            shutil.copy(image_file, tmp_file)
            imagetmp = Image.open(tmp_file)
            imagetmp = imagetmp.rotate(rotate_degree, resample=Image.Resampling.BICUBIC)
            imagetmp.save(tmp_file)
            bio = io.BytesIO()
            imagetmp.save(bio, format='PNG')
            window["image1"].update(data=bio.getvalue(), size=size)
        else:
            continue
        apply_drawing(values, window)
        window["image1"].update(data=bio.getvalue(), size=size)
    '''
    #######################################################################
    # BUTTON APPLY - main event of window
    In this event will be apply the treatment of interferogram image to generate
    the data of gas profile.  
    #######################################################################
    '''
    if event == 'Apply Algorithm':
        # Cleaning plots
        try:
            fig_canvas_agg.get_tk_widget().forget()
        except:
            plt.close('all')
        # Input datas
        h_prof = -1.0
        try:
            begin_x = int(get_value("-BEGIN_X-", values) / scale[0])
            begin_y = int(get_value("-BEGIN_Y-", values) / scale[1])
            end_x = int(get_value("-END_X-", values) / scale[0])
            end_y = int(get_value("-END_Y-", values) / scale[1])
            # get conversion factor in meters/pixel
            factor = float(get_value('-factor-', values)) * 1e-6
            # Manual definition of the filter position
            centerfilter = int(get_value("-centerfilter-", values))
            # sigma value of gaussian function
            sigma = int(get_value('-sigma_gfilter-', values))
            # heat gas constant
            alpha_gas = float(get_value('-polargas-', values)) * 1e-24  # in cm^3
            # Wavelength laser
            lambda0 = float(get_value('-lambda0-', values)) * 1e-9  # in meters
            unc_lambda0 = float(get_value('-unclambda0-', values)) * 1e-9  # in meters
        except:
            sg.popup_error(f"WARNING: Data fields must have numerical values! ")
            continue
        # Colormap distribution
        if values['-cmapcombo-'] == 'Linear':
            colormap_order = 1
        elif values['-cmapcombo-'] == 'Quadratic':
            colormap_order = 2
        elif values['-cmapcombo-'] == 'Cubic':
            colormap_order = 3
        newcmp = func_colormap(colormap_order)

        # Rotate Image
        if rotate_degree != 0:
            originalref = originalref.rotate(total_rot_degree, resample=Image.Resampling.BICUBIC)

        # Input ref. array from ref. image
        array_ref = np.asarray(originalref)

        # Input original files
        phasemaps = []
        gas_dens, gas_abelphasemap, gas_phasemap = [], [], []
        std_phasemap, std_abelmap, std_gas_dens = [],[],[]

        for j in range(0, len(path_files)):
            if rotate_degree != 0:
                originalgas[j] = originalgas[j].rotate(total_rot_degree, resample=Image.Resampling.BICUBIC)
            array_gas = np.asarray(originalgas[j])
            if np.ndim(array_gas) == 3:
                # Slice image with 3 channels:only one channel is used to interferogram treatment
                intref0 = array_ref[:, :, 0]
                intgas0 = array_gas[:, :, 0]
            else:
                intref0 = array_ref[:, :]
                intgas0 = array_gas[:, :]

            # Use whole image or select area of image
            if values['-checkcut-'] == True:
                intref = intref0[begin_y:end_y, begin_x:end_x]
                intgas = intgas0[begin_y:end_y, begin_x:end_x]
            else:
                intref = intref0
                intgas = intgas0

            # Apply Fast Fourier Transform on interferogram data arrays
            fftref = np.fft.fft2(intref)  # ref. interferogram
            fftgas = np.fft.fft2(intgas)  # gas interferogram
            # Defining line or row to apply gaussian filter
            fftmap = np.log(np.abs(fftgas))
            nlmap, nrmap = np.shape(fftmap)
            '''
            # Authomatic definition of the gaussian filter position: 
            this position are defined like the line or column (Vertical or horizontal fringes) with more intensity pixel
            values. This way, this positions are defined using the maximum value of horizontal/vertical pixels sum, 
            depending on fringes orientation.
            Note: Case this filter position is not found, the process is interrupted and the user must select another
            file or another area of interferogram.
            '''
            if centerfilter == 0:
                summap = []
                # Fringe Orientation: HORIZONTAL or VERTICAL
                # sum of array rows (vertical)
                if values['-combofringe-'] == 'vertical':
                    for i in range(0, nrmap):
                        summap.append(np.sum(fftmap[:, i]))
                elif values['-combofringe-'] == 'horizontal':
                    # sum of array lines (horizontal)
                    for i in range(0, nlmap):
                        summap.append(np.sum(fftmap[i,]))

                # Defining point of gaussian filter application using max value of horizontal/vertical pixels sum
                filterpoints, _ = find_peaks(summap, height=0.9 * np.max(summap))
                # Range of gaussian filter
                filterspoints_widths = (peak_widths(summap, filterpoints, rel_height=0.5)[0])

                if len(filterpoints) == 0:
                    sg.popup(f"WARNING: Unable to apply the Fast Fourier Transform to the selected image!")
                    continue
                if values['-combofringe-'] == 'horizontal':
                    # filter range is equal to FWHM of signal of summaps
                    if filterpoints[0] <= 5:
                        centerfilter = filterpoints[1]
                        f_range = int(filterspoints_widths[1])
                    else:
                        centerfilter = filterpoints[0]
                        f_range = int(filterspoints_widths[0])
                elif values['-combofringe-'] == 'vertical':

                    if filterpoints[len(filterpoints) - 1] >= nrmap - 5:
                        centerfilter = filterpoints[len(filterpoints) - 2]
                        f_range = int(filterspoints_widths[len(filterpoints) - 2])
                    else:
                        centerfilter = filterpoints[len(filterpoints) - 1]
                        f_range = int(filterspoints_widths[len(filterpoints) - 1])

            window['-centerfilter-'].update(str(centerfilter))

            # Creating filter array from null array
            gfilter = np.zeros(np.shape(fftgas))

            # Creating Filter for Horizontal/vertical fringes orientation
            if values['-combofringe-'] == 'horizontal':
                gfilter[centerfilter - f_range:centerfilter + f_range] = np.ones(np.shape(
                    fftgas[centerfilter - f_range:centerfilter + f_range]))
                # Applying gaussian filter at selected filter position
                gfilter = gaussian_filter(gfilter, sigma=2 * f_range)
            elif values['-combofringe-'] == 'vertical':
                gfilter[:, centerfilter - f_range:centerfilter + f_range] = np.ones(np.shape(
                    fftgas[:, centerfilter - f_range:centerfilter + f_range]))
                # Applying gaussian filter at selected filter position
                gfilter = gaussian_filter(gfilter, sigma=2 * f_range)

            # Applying Inverse FFT in resultant array obtained after use of the gaussian filter on FFT arrays
            ifftref = np.fft.ifft2(gfilter * fftref)
            ifftgas = np.fft.ifft2(gfilter * fftgas)

            # Creating Phase Maps arrays by subtracting the arguments of IFFT arrays
            phasemaps = (np.angle(ifftgas) - np.angle(ifftref))
            # Unwrap phase:
            uwphasemap = unwrap_phase(phasemaps)
            # Range for scan is 5% of total dimension of matrix
            frgs_shifts = np.ones(np.shape(intref))
            frgs_widths = np.ones(np.shape(intref))
            frgs_shifts, frgs_widths = [], []
            lim_l = int(0.05 * nlmap) + 1
            lim_r = int(0.05 * nrmap) + 1
            if values['-combofringe-'] == 'vertical':
                # Intensity distribution
                dist1 = intensity_dist(intref, 0)
                dist2 = intensity_dist(intgas, 0)
                for l in range(0, lim_l):
                    # Scanning fringes pattern to define min. shift
                    frgs_ref = (intref0[l, begin_x:end_x])
                    frgs_plasma = (intgas0[l, begin_x:end_x])
                    frgs_shifts_i, frgs_widths_i = fringes_info(frgs_plasma, frgs_ref, intref)
                    frgs_shifts.append(frgs_shifts_i)
                    frgs_widths.append(frgs_widths_i)

            if values['-combofringe-'] == 'horizontal':
                # Intensity distribution
                dist1 = intensity_dist(intref, 1)
                dist2 = intensity_dist(intgas, 1)
                # fringes shift
                for l in range(0, lim_r):
                    frgs_ref = np.transpose(intref0[begin_y:end_y, l])
                    frgs_plasma = np.transpose(intgas0[begin_y:end_y, l])
                    frgs_shifts_i, frgs_widths_i = fringes_info(frgs_plasma, frgs_ref, np.transpose(intref))
                    frgs_shifts.append(np.transpose(frgs_shifts))
                    frgs_widths.append(np.transpose(frgs_widths))

            frgs_shifts = mean_maps(frgs_shifts)
            frgs_widths = mean_maps(frgs_widths)
            frgs_shifts = gaussian_filter(frgs_shifts, sigma=int(np.mean(frgs_widths) / 2))
            frgs_widths = gaussian_filter(frgs_widths, sigma=int(np.mean(frgs_widths) / 2))
            std_phasemap_i = ((np.pi * frgs_shifts / (2 * frgs_widths)) * \
                              np.sqrt((np.mean(dist1) * (dist1 + dist2)) / (2 * dist1 * dist2)))
            std_phasemap.append(std_phasemap_i)

            '''
            ################################################################################
            Applying Inverse Abel Transform (IAT):
            The IAT is applied using PyAbel algorithm and to apply its library correctly is necessary
            to define a axis symmetric in image (Horizontal or Vertical). In gas profile case the axissymmetric is 
            defined by more intensity pixel range. So, the image is cut according axissymetric.
            The right side of image is used like standard to use IAT.
            NOTE: the Abel transform is always performed around the vertical axis, so when the image have horizontal
            axissymmetry the matrix must be transposed.

            '''
            # Transpose Matrix for Horizontal Axissmetry
            if values['-comboaxisymm-'] == 'horizontal':
                phasemap_corr = np.transpose(phasemap_corr)
            # Apply gaussian filter to define the region with more intensity pixel value
            phasemap_corr = gaussian_filter(uwphasemap, sigma=sigma)
            phasemap_corr = phasemap_corr - np.ones(np.shape(phasemap_corr))*np.min(phasemap_corr)
            nlines, nrows = np.shape(phasemap_corr)

            # Define region with more intensity pixel - position x and y
            cline, crow = np.where(phasemap_corr >= phasemap_corr.max() * 0.95)
            cy, cx = int(np.median(cline)), int(np.median(crow))
            vert_lim = int(2*cx+1)
            # If the region not found, set symmetric point like half image
            if math.isnan(cx) == True:
                cx = nrow // 2

            phasemap_symm = phasemap_corr[:, 0:vert_lim]
            gas_diameter = int(2 * cy) # diameter of the cylinder used in IAT

            try:
                # Applying inverse Abel Transform
                phase_abel = abel.Transform((phasemap_symm), symmetry_axis=0, direction='inverse',
                                              method='onion_peeling').transform
            except:
                phase_abel = np.zeros(np.shape(phasemap_symm))
                sg.popup_error(f"WARNING: Unable to apply the Abel transform to the selected image! ")

            '''
            ############################################################################################
            Calculating std from Abel Transform:
            The std is calculated using deviation of mormalized phasemap and normalized IAT phasemap  
            '''
            rangeh0, rangev0 = np.shape(phase_abel)
            norm_phasemap = np.zeros(np.shape(phase_abel))
            for k in range(0, rangeh0):
                norm_phasemap[k] = phasemap_symm[k] * np.max(abs(phase_abel[k])) / np.max(abs(phasemap_symm[k]))
            std_abelmap0 = np.sqrt(np.square(phase_abel-norm_phasemap))

            if values['-comboaxisymm-'] == 'horizontal':
                vert_lim = nlines
                phasemap_corr = np.transpose(phasemap_corr)
                phasemap_symm = np.transpose(phasemap_symm)
                phase_abel = np.transpose(phase_abel)
                std_abelmap0 = np.transpose(std_abelmap0)
            gas_phasemap.append(phasemap_corr)
            '''
            ########################################################################################
            Calculating refraction index and gas density from IAT phasemap:
            This step is subdivided in 2 modes of refrac. index and gas dens. calculation.
            1 - Symmetric Images:
                In this case refrac. index and gas density are calculated directly by IAT phasemap,
                exploring the symmetry of image

            2 - For No Symmetric Images
                In this case the refrac. index is calculated directly by phasemap (without IAT) and
                gas dens. are normalized using gas density values obtained using IAT for each height
                of the gas density profile.

            NOTE: The gas density is calculated using C-M relation.        
            '''
            # Calculating index refraction from IAT of phasemap
            n_index0 = (1 + (phase_abel * lambda0) / (2 * np.pi*factor))
            # Cutting border of images due the computational artefacts generated by IAT and problems with no symmetric images
            n_index = n_index0[:, int(0.05 * vert_lim): int(0.95 * vert_lim)]
            # Calculating gas density from C-M relation
            gas_dens_i = (3 * (np.square(n_index) - np.ones(np.shape(n_index)))) / \
                             (4 * np.pi * alpha_gas * (np.square(n_index) + 2 * np.ones(np.shape(n_index))))

            #new matrix size for plot
            rangeh, rangev = np.shape(gas_dens_i)

            std_abelmap_i = std_abelmap0[:, int(0.05 * vert_lim): int(0.95 * vert_lim)]
            std_abelmap.append(std_abelmap_i)

            gas_abelphasemap.append(phase_abel)
            gas_dens.append(gas_dens_i)

            '''
            CALCULATION TOTAL STANDARD DEVIATION FROM:
            1. Measurement of interferogram
            2. Inverse Abel Transform
            3. Laser wavelength
            '''
            dN_n = ((9 / (2 * np.pi * alpha_gas)) * n_index) / np.square(
                np.square(n_index) + 2 * np.ones(np.shape(n_index)))
            # Contribution 1: measurement interferogram
            std_phase1 = np.square(std_phasemap_i[:, int(0.05 * vert_lim): int(0.95 * vert_lim)]\
                                   /(factor*gas_diameter)) #rad/metro
            # Contribution 2: Abel transformation accuracy
            std_phase2 = np.square(std_abelmap_i/factor) #rad/metro
            std_phase = np.sqrt(std_phase2+std_phase1) #rad/metro
            dn_phase=(lambda0/(2*np.pi)) #m

            # Contribution 3: laser wavelength
            dn_lambda = (phase_abel[:, int(0.05 * vert_lim): int(0.95 * vert_lim)]\
                                  /(2*np.pi*factor))#rad/m

            std_gas_dens_i = np.sqrt(np.square(dN_n)*(np.square(dn_phase*std_phase) + \
                                                      np.square(dn_lambda*unc_lambda0)))*1e-6 # cm^-3
            std_gas_dens.append(std_gas_dens_i)

        #BUILDING MATRIX RESULTS FOR:
            # Many files
        if len(gas_phasemap)>1:
            # PHASEMAP
            gas_phasemap_mean = mean_maps(gas_phasemap)
            std_phasemap_mean = np.sqrt(np.square(mean_maps(std_phasemap))+\
                                                  np.square(std_maps(gas_phasemap,gas_phasemap_mean)))
            #INV. ABEL TRANSF. MAP
            gas_abelmap_mean = mean_maps(gas_abelphasemap)
            std_abelmap_mean = np.sqrt(np.square(mean_maps(std_abelmap)) + \
                                        np.square(std_maps(gas_abelphasemap, std_abelmap_mean)))
            # PLASMA DENSITY
            gas_dens_mean = mean_maps(gas_dens)
            std_dens_mean = np.sqrt(np.square(mean_maps(std_gas_dens)) + \
                                       np.square(std_maps(gas_dens, std_dens_mean)))
        else:
            #PHASEMAP
            gas_phasemap_mean = (gas_phasemap[0])
            std_phasemap_mean =(std_phasemap[0])
            #INV. ABEL TRANSF. MAP
            gas_abelmap_mean = (gas_abelphasemap[0])
            std_abelmap_mean = (std_abelmap[0])
            #PLASMA DENSITY
            gas_dens_mean = (gas_dens[0])
            std_dens_mean = (std_gas_dens[0])

        '''
        BUILDING 2D AND 1D PLOTS
        '''
        #Ajust slider for horizontal/vertical
        if values['-comboaxisymm-'] == 'vertical':
            window['sliderh'].update(range=(0, rangev - 1))
        elif values['-comboaxisymm-'] == 'horizontal':
            window['sliderh'].update(range=(0, rangeh - 1))

        # Plots are building from user select
        if values['fftradio'] == True:  # Plot FFT map result
            matrix_plot = fftmap
        elif values['filterradio'] == True:  # Plot gaussian filter map
            matrix_plot = gfilter
        elif values['phaseradio'] == True:  # Plot phase map result
            if values['-checkstd-'] == False:
                matrix_plot = gas_phasemap_mean
            else:
                matrix_plot = std_phasemap_mean
        elif values['abelradio'] == True:  # Plot gas density profile from IAT
            if values['-checkstd-'] == False:
                matrix_plot = gas_abelmap_mean
            else:
                matrix_plot = std_abelmap_mean
        elif values['densradio'] == True:  # Plot gas density profile
            if values['-checkstd-'] == False:
                matrix_plot = gas_dens_mean# - np.ones(np.shape(plasma_dens_mean)) * np.min(plasma_dens_mean)
            else:
                matrix_plot = std_dens_mean

        # colormap distribution definition
        if values['-cmapcombo-'] == 'Linear':
            colormap_order = 1
        elif values['-cmapcombo-'] == 'Quadratic':
            colormap_order = 2
        elif values['-cmapcombo-'] == 'Cubic':
            colormap_order = 3

        # Creating plot figure parameters
        fig, ax1 = plt.subplots(figsize=(4.9, 4))

        if values['filterradio'] == True:
            abel_plot = ax1.imshow(matrix_plot, cmap='gray')

        elif values['fftradio'] == True:
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
            if values['-checkstd-'] == False:
                cb1.set_label(label='$Accumulatad\hspace{.5}Phase\hspace{.5} (rad)$', size=12, weight='bold')
            else:
                cb1.set_label(label='$Standard\hspace{.5}dev.\hspace{.5}Acc.\hspace{.5}Phase\hspace{.5} (rad)$',\
                              size=12, weight='bold')

        elif values['densradio'] == True:
            divider = make_axes_locatable(ax1)
            extentplot = np.shape(matrix_plot)
            x_max = extentplot[1] * factor * 1e6
            y_max = extentplot[0] * factor * 1e6
            cax = divider.append_axes("right", size="5%", pad=0.05)
            abel_plot = ax1.imshow(matrix_plot, extent=[0, x_max, 0, y_max], cmap=newcmp)
            cb1 = fig.colorbar(abel_plot, cax=cax)
            ax1.set_xlabel('$x\hspace{.5}(\mu m)$', fontsize=12)
            ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
            if values['-checkstd-'] == False:
                cb1.set_label(label='$Gas\hspace{.5}Density\hspace{.5} (cm^{-3})$', \
                              size=12, weight='bold')
            else:
                cb1.set_label(label='$Standard\hspace{.5}dev.\hspace{.5}Gas\hspace{.5}Density\hspace{.5} (cm^{-3})$',\
                              size=12, weight='bold')

        else:
            divider = make_axes_locatable(ax1)
            extentplot = np.shape(matrix_plot)
            x_max = extentplot[1] * factor * 1e6
            y_max = extentplot[0] * factor * 1e6
            cax = divider.append_axes("right", size="5%", pad=0.05)

            abel_plot = ax1.imshow(matrix_plot, extent=[0, x_max, 0, y_max], cmap=newcmp)
            cb1 = fig.colorbar(abel_plot, cax=cax)
            cb1.set_label(label='$Gas\hspace{.5}Density\hspace{.5} (cm^{-3})$', size=12, weight='bold')
            ax1.set_xlabel('$x\hspace{.5}(\mu m)$', fontsize=12)
            ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
            if values['-checkstd-'] == False:
                cb1.set_label(label='$Phase\hspace{.5}Map\hspace{.5} (rad)$', size=12, weight='bold')
            else:
                cb1.set_label(label='$Standard\hspace{.5}dev.\hspace{.5}Phase\hspace{.5}Map\hspace{.5} (rad)$',\
                              size=12, weight='bold')

        fig.tight_layout(pad=2)
        fig_canvas_agg = draw_figure(window['canvasabel'].TKCanvas, fig)

        visible_f1d = False
        # Enable/Disable specific buttons and frames for 2D analysis
        window['Save Data'].update(disabled=False)
        window['Save Plot'].update(disabled=False)
        window['frame1d'].update(visible=False)
        window['Dens. Profile 2D'].update(disabled=False)
        window['Dens. Profile 1D'].update(disabled=False)
    #########################################################################
    # BUTTON DENS.PROFILE 2D
    if event == 'Dens. Profile 2D':
        # set height position
        h_prof = -1.0
        # Plots are building from user select
        if values['fftradio'] == True:  # Plot FFT map result
            matrix_plot = fftmap
        elif values['filterradio'] == True:  # Plot gaussian filter map
            matrix_plot = gfilter
        elif values['phaseradio'] == True:  # Plot phase map result
            if values['-checkstd-'] == False:
                matrix_plot = gas_phasemap_mean
            else:
                matrix_plot = std_phasemap_mean
        elif values['abelradio'] == True:  # Plot gas density profile from IAT
            if values['-checkstd-'] == False:
                matrix_plot = gas_abelmap_mean
            else:
                matrix_plot = std_abelmap_mean
        elif values['densradio'] == True:  # Plot gas density profile
            if values['-checkstd-'] == False:
                matrix_plot = gas_dens_mean - np.ones(np.shape(gas_dens_mean)) * np.min(gas_dens_mean)
            else:
                matrix_plot = std_dens_mean

        # colormap distribution definition
        if values['-cmapcombo-'] == 'Linear':
            colormap_order = 1
        elif values['-cmapcombo-'] == 'Quadratic':
            colormap_order = 2
        elif values['-cmapcombo-'] == 'Cubic':
            colormap_order = 3
        newcmp = func_colormap(colormap_order)

        try:
            # clearing figures and plots
            fig_canvas_agg.get_tk_widget().forget()
            plt.close('all')

            # Instead of plt.show
            fig, ax1 = plt.subplots(figsize=(4.9, 4))

            if values['filterradio'] == True:
                abel_plot = ax1.imshow(matrix_plot, cmap='gray')

            elif values['fftradio'] == True:
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
                if values['-checkstd-'] == False:
                    cb1.set_label(label='$Accumulatad\hspace{.5}Phase\hspace{.5} (rad)$', size=12, weight='bold')
                else:
                    cb1.set_label(label='$Standard\hspace{.5}dev.\hspace{.5}Acc.\hspace{.5}Phase\hspace{.5} (rad)$', \
                                  size=12, weight='bold')

            elif values['densradio'] == True:
                divider = make_axes_locatable(ax1)
                extentplot = np.shape(matrix_plot)
                x_max = extentplot[1] * factor * 1e6
                y_max = extentplot[0] * factor * 1e6
                cax = divider.append_axes("right", size="5%", pad=0.05)
                abel_plot = ax1.imshow(matrix_plot, extent=[0, x_max, 0, y_max], cmap=newcmp)
                cb1 = fig.colorbar(abel_plot, cax=cax)
                ax1.set_xlabel('$x\hspace{.5}(\mu m)$', fontsize=12)
                ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
                if values['-checkstd-'] == False:
                    cb1.set_label(label='$Plasma\hspace{.5}Density\hspace{.5} (cm^{-3})$', \
                                  size=12, weight='bold')
                else:
                    cb1.set_label(label='$Standard\hspace{.5}dev.\hspace{.5}Plasma\hspace{.5}Density\hspace{.5} (cm^{-3})$', \
                                  size=12, weight='bold')

            else:
                divider = make_axes_locatable(ax1)
                extentplot = np.shape(matrix_plot)
                x_max = extentplot[1] * factor * 1e6
                y_max = extentplot[0] * factor * 1e6
                cax = divider.append_axes("right", size="5%", pad=0.05)

                abel_plot = ax1.imshow(matrix_plot, extent=[0, x_max, 0, y_max], cmap=newcmp)
                cb1 = fig.colorbar(abel_plot, cax=cax)
                cb1.set_label(label='$Plasma\hspace{.5}Density\hspace{.5} (cm^{-3})$', size=12, weight='bold')
                ax1.set_xlabel('$x\hspace{.5}(\mu m)$', fontsize=12)
                ax1.set_ylabel('$y\hspace{.5}(\mu m)$', fontsize=12)
                if values['-checkstd-'] == False:
                    cb1.set_label(label='$Phase\hspace{.5}Map\hspace{.5} (rad)$', size=12, weight='bold')
                else:
                    cb1.set_label(label='$Standard\hspace{.5}dev.\hspace{.5}Phase\hspace{.5}Map\hspace{.5} (rad)$', \
                                  size=12, weight='bold')

            fig.tight_layout(pad=2)
            fig_canvas_agg = draw_figure(window['canvasabel'].TKCanvas, fig)

            visible_f1d = False
            window['frame1d'].update(visible=False)

        except:
            continue

    #########################################################################
    # BUTTON DENS.PROFILE 1D AND SLIDER POSITION
    if (event == 'Dens. Profile 1D') or (event == 'sliderh'):
        # set height position
        h_prof = -1.0
        # Plots are building from user select
        if values['fftradio'] == True:  # Plot FFT map result
            matrix_plot = fftmap
            matrix_plot_std = np.zeros(np.shape(matrix_plot))
        elif values['filterradio'] == True:  # Plot gaussian filter map
            matrix_plot = gfilter
            matrix_plot_std = np.zeros(np.shape(matrix_plot))
        elif values['phaseradio'] == True:  # Plot phase map result
            if values['-checkstd-'] == False:
                matrix_plot = gas_phasemap_mean
                matrix_plot_std = std_phasemap_mean
        elif values['abelradio'] == True:  # Plot gas density profile from IAT
            matrix_plot = gas_abelmap_mean
            matrix_plot_std = std_abelmap_mean
        elif values['densradio'] == True:  # Plot gas density profile
            matrix_plot = gas_dens_mean
            matrix_plot_std = std_dens_mean

        try:
            # clearing figures and plots
            fig_canvas_agg.get_tk_widget().forget()
            plt.close('all')

            # create r axis according to symmetry in micrometers (µm)
            rangeh, rangev = np.shape(matrix_plot)
            if values['-comboaxisymm-'] == 'vertical':
                raxis = np.arange(-rangev / 2, rangev / 2, 1)
                raxis_um = raxis * factor * 1e6  # um
                # set origin position (exit nozzle position) and slider position
                h_prof = 0
                pos_0 = rangeh - 1
                pos = pos_0 - int(values['sliderh'])
                # convert vertical array positions to height positions in µm
                h_prof = int(values['sliderh']) * factor * 1e6
                array_plot = matrix_plot[pos]
                array_std = matrix_plot_std[pos]

            elif values['-comboaxisymm-'] == 'horizontal':
                window['sliderh'].update(range=(0, rangeh - 1))
                raxis = np.arange(-rangeh / 2, rangeh / 2, 1)
                raxis_um = raxis * factor * 1e6  # um
                # set origin position (exit nozzle position) and slider position
                h_prof = 0
                pos = int(values['sliderh'])
                # convert vertical array positions to height positions in µm
                h_prof = int(values['sliderh']) * factor * 1e6
                array_plot = matrix_plot[:, pos]
                array_std = matrix_plot_std[:, pos]

            # Creating plot parameters
            fig, ax1 = plt.subplots(figsize=(4.9, 4))
            ax1.plot(raxis_um, array_plot, label='$%d \hspace{.5}\mu m$' % h_prof, lw=2, color="blue")
            if values['-checkstd-'] == True:
                # ax1.fill_between(raxis_um, array_plot-array_std,array_plot+array_std,label='$\sigma_{dens.}$',
                # alpha=0.2, color="blue" )
                ax1.set_ylim(0., np.max(array_std + array_plot) * 1.05)
                ax1.errorbar(raxis_um, array_plot, yerr=array_std, label='$\sigma_{dens.}$', alpha=0.2,
                             color="blue")

            # Including new 1D density profile for another height from origin height position
            if values['-checkpos1-'] == True:
                h_prof1 = int(get_value('-pos1-', values))
                pos1 = int(h_prof1 / (factor * 1e6))
                if values['-comboaxisymm-'] == 'vertical':
                    ax1.plot(raxis_um, matrix_plot[pos1], label='$%d \hspace{.5}\mu m$' % h_prof1, lw=1,
                             color="red")
                elif values['-comboaxisymm-'] == 'horizontal':
                    ax1.plot(raxis_um, matrix_plot[:, pos1], label='$%d \hspace{.5}\mu m$' % h_prof1, lw=1,
                             color="red")
            if values['-checkpos2-'] == True:
                h_prof2 = int(get_value('-pos2-', values))
                pos2 = int(h_prof2 / (factor * 1e6))
                if values['-comboaxisymm-'] == 'vertical':
                    ax1.plot(raxis_um, matrix_plot[pos2], label='$%d \hspace{.5}\mu m$' % h_prof2, lw=1,
                             color="red")
                elif values['-comboaxisymm-'] == 'horizontal':
                    ax1.plot(raxis_um, matrix_plot[:, pos2], label='$%d \hspace{.5}\mu m$' % h_prof2, lw=1,
                             color="red")
            if values['-checkpos3-'] == True:
                h_prof3 = int(get_value('-pos3-', values))
                pos3 = int(h_prof3 / (factor * 1e6))
                if values['-comboaxisymm-'] == 'vertical':
                    ax1.plot(raxis_um, matrix_plot[pos3], label='$%d \hspace{.5}\mu m$' % h_prof3, lw=1,
                             color="red")
                elif values['-comboaxisymm-'] == 'horizontal':
                    ax1.plot(raxis_um, matrix_plot[:, pos3], label='$%d \hspace{.5}\mu m$' % h_prof3, lw=1,
                             color="red")

            ax1.set_xlabel('$r\hspace{.5}(\mu m)$', fontsize=12)

            if values['densradio'] == True:
                ax1.set_ylabel('$Plasma\hspace{.5}Density\hspace{.5} (cm^{-3})$', fontsize=12)

            elif values['abelradio'] == True:
                ax1.set_ylabel('$Phase\hspace{.5} (rad)$', fontsize=12)

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
        save_filename_data = sg.popup_get_file('File', file_types=[("DAT (*.dat)", "*.dat"), ("TXT (*.txt)", "*.txt")],
                                               save_as=True, no_window=True)

        if save_filename_data:
            file_data = open(save_filename_data, 'a')
            file_data.seek(0)  # sets  point at the beginning of the file
            file_data.truncate()

            if visible_f1d==True:
                rangeh, rangev = np.shape(matrix_plot)
                raxis = np.arange(-rangev / 2, rangev / 2, 1)
                raxis_um = raxis * factor * 1e6  # um
                # save data plot
                file_data.write('\nx[µm] \t f(x)\n')
                # Plot1D
                if h_prof >= 0.0:
                    file_data.write('\t (%.0f µm)' % h_prof)
                    list_data = np.vstack((raxis_um, array_plot))
                    if values['-checkstd-'] == True:
                        file_data.write('\t (error)')
                        list_data = np.vstack((raxis_um, array_std))
                    # verify additional height positions 1, 2 and 3
                    if values['-checkpos1-'] == True:
                        if values['-comboaxisymm-'] == 'vertical':
                            list_data = np.vstack((list_data, matrix_plot[pos1]))
                        else:
                            list_data = np.vstack((list_data, matrix_plot[:, pos1]))
                        file_data.write('\t(%.0f µm)' % h_prof1)

                    if values['-checkpos2-'] == True:
                        if values['-comboaxisymm-'] == 'vertical':
                            list_data = np.vstack((list_data, matrix_plot[pos2]))
                        else:
                            list_data = np.vstack((list_data, matrix_plot[:, pos2]))
                        file_data.write('\t(%.0f µm)' % h_prof2)

                    if values['-checkpos3-'] == True:
                        if values['-comboaxisymm-'] == 'vertical':
                            list_data = np.vstack((list_data, matrix_plot[pos3]))
                        else:
                            list_data = np.vstack((list_data, matrix_plot[:, pos3]))
                        file_data.write('\t(%.0f µm)' % h_prof3)

                    file_data.write('\n')
                    list_str = (np.transpose(list_data))
                    np.savetxt(file_data, list_str, fmt='%.2e', delimiter='\t')

                elif h_prof < 0.0:
                    # Write density in H=0

                    list_data = np.vstack((raxis_um, matrix_plot))
                    file_data.write('\nx[µm] \t\t f(x)\n')
                    list_str = (np.transpose(list_data))
                    np.savetxt(file_data, list_str, fmt='%.2e')
                sg.popup(f"Saved: {save_filename_data}")
                file_data.close()
            else:
                if values['-checkstd-'] == False:
                    np.savetxt(file_data, matrix_plot, fmt='%.3e')
                else:
                    np.savetxt(file_data, matrix_plot_std, fmt='%.3e')
                sg.popup(f"Saved: {save_filename_data}")
                file_data.close()
        else:
            continue

window.close()