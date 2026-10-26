#!/usr/bin/python3

# if config file doesn't exist, make a template. 
def make_config_template():
    config_file.write(
        "Compounds: []",
        "#Example, Compounds: [a,b,c]",
        "Molecular weights: []",
        "#Example, Molecular weights: [10,42,67]",
        "#The molecular weight can be a range, use format [default,floor,ceiling]",
        "Stock solution volume: []",
        "#Example, in units of ml Stock solution volume: [25,25,25]"
        "Stock solution compound mass as measured: []",
        "#Example, [1,2,[1,1,1]] means 1g of compound a for series 1, 2 g of compound b for series 2, 1 g of compound a,b,c each for series 3.",
        "Solution origins: []",
        "#Need one List for each compound. S means stock.",
        "#[S,2,1] means solution is made from 2 ml of stock and 1 ml of water.",
        "#[-1,2,1] means solution is made from 2 ml of the previous solution, 1 ml of water",
        "#[4,0.2,2.8] means the solution is made from 0,2 ml of the 4th entry and 2.8 ml of water",
        "#[[S,None,S],[0.5,0.5],2] means the solution is made from the stock solution of the first and third compound at 0,5 ml each, with 2 ml of added water.",
        "Primary peaks: [None]",
        "Peak criteria: [None]",
        "#If primary peak is given but not peak criteria, the script only looks at that wavelength for the corresponding compound.",
        "#Peak criteria is the fraction of the strongest peak. If primary peaks is not given but the criteria is given, it will use all peaks above the criteria.",
        "#If neither is given, any peak that's at 0.5 of the strongest peak or higher will be used."
        )
        
# read config file
def read_config():
    #  calcu_concentrations is a list. 
    #Each element in the list corresponds to an element in experiment. 
    #Each element in the list contains the calculated concentration for each solution.
    return calc_concentrations

# for finding peaks when primary peaks is not given.
def find_peaks(peak_criteria):
    #for now, don't implement this because we know what peaks we want to use.
    return None

# return the absorbance at the specified wavelength.
def read_absorbance(primary_peaks):
    return absorbance

# determine linear range
def determine_linear_range(spectrum):
    return linear_range_floor,linear_range_ceiling

# create calibration curve for a single compound series
def single_calibration_curve(absorbances,concentrations):
    return a,b, rsq # y = ax + b

def mixture_calibration_curve(mix_absorbances,mix_concentrations):
    return slopes,intersects,rsqs # y1 = ax1 + b, y2 = cx2 + d, assuming mixture absorbace y = y1 + y2

def deconvolute(spectrum,concentrations,single_calibration_curves):
    return compositions,rsq

# determine if we should use deconvolution method or we should use calibration curve from mixture data, or if we should just interpolate using nearest data points
def determine_analysis_method():
    return use_deconvolution,use_mix_calibration_curve,use_nearest_interpolation

# just to look at a scatter plot of concentration and absorbance 
def data_visualization(compond,series):
    # Make a figure for each of the compound for each solution series that only contain one compound

    # Make a figure for all coumpounds in a mixture 

    # Make a figure for each compound showing all data, but show data points from mixture series using a different color.

    return None

# analyze the analate.
def analyze(analate_spectrum):

    return concentrations


