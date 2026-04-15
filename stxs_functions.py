import importlib.util
spec = importlib.util.find_spec("vbfprocessor")
print(spec.origin)

import sys
sys.path.append('/uscms/home/azhou/nobackup/smeft/analysis/hbb-coffea/')
print(sys.path)
import re
import json
from scipy.optimize import curve_fit
import hist
import vbfprocessor
import awkward as ak
import uproot
import os
import time
import matplotlib.pyplot as plt
import numpy as np

import coffea
print(coffea.__version__)
from coffea import util, processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, PFNanoAODSchema
from coffea.processor import Runner, FuturesExecutor, IterativeExecutor
import boostedhiggs.ewk_higgs_correction as ewk
importlib.reload(ewk)
print("vbfprocessor path:", vbfprocessor.__file__)
importlib.reload(vbfprocessor)
#############################################

def get_wc_order(reweight_card, dictionary_path='dictionary.json'):
    with open(dictionary_path) as f:
        d = json.load(f)
    order = []
    seen = set()
    with open(reweight_card) as f:
        for line in f:
            if line.startswith('set '):
                _, model, param_index, _ = line.split()
                name = d[model][param_index]
                if name not in seen:
                    order.append(name)
                    seen.add(name)
    return order
#print(get_wc_order('/uscms/home/azhou/nobackup/smeft/cmseft/generation/genproductions/bin/MadGraph5_aMCatNLO/cards/VBF_SMEFTsim_topU3l_NP1/VBF_SMEFTsim_topU3l_NP1_reweight_card.dat'))

def generate_tuples(n):
    values = [0.0, 0.5, 1.0]
    tuples = [()]
    for _ in range(n):
        tuples = [t + (v,) for t in tuples for v in values]
    return tuples

def wc_map_dict(unsorted_wc_list):
    
    wc_order = get_wc_order('/uscms/home/azhou/nobackup/smeft/cmseft/generation/genproductions/bin/MadGraph5_aMCatNLO/cards/VBF_SMEFTsim_topU3l_NP1/VBF_SMEFTsim_topU3l_NP1_reweight_card.dat')
    invalid_operators = [x for x in unsorted_wc_list if x not in wc_order]
    if invalid_operators:
        raise ValueError(f"WCs not found in wc_order built on reweight card (see path): {invalid_operators}.")
    
    wc_list = sorted(unsorted_wc_list, key=lambda x: wc_order.index(x))
    
    point_list = generate_tuples(len(wc_list))
    wc_map = {}
    for point in point_list:
        non_sm_points = [(wc_list[i], point[i]) for i in range(len(point)) if point[i] != 0.0]
        if not non_sm_points:
            name = "SM"
        else:
            name = ",".join(f"{wc}={val}" for wc, val in non_sm_points)
        wc_map[point] = name
    
    return wc_map

def get_bin_yields(hslice):
    bin1_codes = [221, 222]
    bin2_codes = [223, 224]
    
    axis = hslice.axes["htxs_stage2"]
    values = hslice.values()
    codes = [axis.bin(i) for i in range(axis.size)]
    code_to_value = {code: values[i] for i, code in enumerate(codes)}
    bin1 = sum(code_to_value.get(c, 0.0) for c in bin1_codes)
    bin2 = sum(code_to_value.get(c, 0.0) for c in bin2_codes)
    return bin1, bin2

def quad_2D(xdata, a1, a2, b1, b2, c, d):
    x1, x2 = xdata
    return a1*x1*x1 + a2*x1 + b1*x2*x2 + b2*x2 + c*x1*x2 + d

def quad_1D(xdata, a1, a2, a3):
    x1 = xdata
    return a1*x1*x1 + a2*x1 + a3

def stxs_reweight_function(coffea_name, operator_list, verbose = False):
    
    bin1_codes = [221, 222]
    bin2_codes = [223, 224]

    bin1_yield_list = {}
    bin2_yield_list = {}

    all_h = util.load(f"coffea/{coffea_name}.coffea")["htxs"]
    
    wc_mapping = wc_map_dict(operator_list)
    if verbose == True:
        print('wc_mapping:', wc_mapping)

    for p, wc_label in wc_mapping.items():
        if wc_label not in list(all_h.axes["wc"]):
            continue
        hslice = all_h[{"wc": wc_label}]
        bin1, bin2 = get_bin_yields(hslice)
        bin1_yield_list[p] = bin1
        bin2_yield_list[p] = bin2

    if verbose == True:
        print('bin1 yield dict: ', bin1_yield_list)
        print('bin2 yield dict: ',bin2_yield_list)

    points = wc_mapping.keys()
    points = [point for point in points if point in bin1_yield_list.keys() and point in bin2_yield_list]
    
    y_bin1 = np.array([bin1_yield_list[point] for point in points])
    y_bin2 = np.array([bin2_yield_list[point] for point in points])


    if len(operator_list) == 2:
        if len(points) < 6:
            raise ValueError("Need at least 6 points for 2D quadratic fit")
        op1_vals = np.array([point[0] for point in points])
        op2_vals = np.array([point[1] for point in points])
        
        fit_dim = int(((len(operator_list) + 1) * (len(operator_list) + 2))/2)
        if verbose == True:
            print('fit_dim:', fit_dim)
        p0 = list(np.ones(fit_dim, dtype=float)) #initial guess
        if verbose == True:
            print(p0)
        
        coeff_bin1, cov_bin1 = curve_fit(quad_2D, (op1_vals, op2_vals), y_bin1, p0=p0)
        coeff_bin2, cov_bin2 = curve_fit(quad_2D, (op1_vals, op2_vals), y_bin2, p0=p0)

    elif len(operator_list) == 1:
        if len(points) < 3:
            raise ValueError("Need at least 3 points for 1D quadratic fit")
            
        op1_vals = np.array([point[0] for point in points])
        
        fit_dim = int(((len(operator_list) + 1) * (len(operator_list) + 2))/2)
        if verbose == True:
            print('fit_dim:', fit_dim)
        p0 = list(np.ones(fit_dim, dtype=float)) #initial guess
        if verbose == True:
            print(p0)
        
        coeff_bin1, cov_bin1 = curve_fit(quad_1D, (op1_vals), y_bin1, p0=p0)
        coeff_bin2, cov_bin2 = curve_fit(quad_1D, (op1_vals), y_bin2, p0=p0)

    if verbose == True:
        print('coeff_bin1: ', coeff_bin1)
        print('coeff_bin2: ', coeff_bin2)

    return coeff_bin1, coeff_bin2

def chisq_stxs(calculated_xsec1, calculated_xsec2):
    cms_obs = np.array([240., 120., 200., 190., 68., 61.])
    bin1_obs = cms_obs[0]
    bin2_obs = cms_obs[1]

    sig1_up = cms_obs[2]
    sig1_down = cms_obs[3]
    sig2_up = cms_obs[4]
    sig2_down = cms_obs[5]
    
    predicted = np.array([calculated_xsec1,calculated_xsec2])
    
    diff = predicted - cms_obs[0:2]
    bin1_diff = diff[0]
    bin2_diff = diff[1]
    
    sigma1 = sig1_up if bin1_diff >=0 else sig1_down
    sigma2 = sig2_up if bin2_diff >= 0 else sig2_down

    bin1_chisq = ((calculated_xsec1 - bin1_obs) / sigma1) **2
    bin2_chisq = ((calculated_xsec2 - bin2_obs) / sigma2) **2
    total_chisq = bin1_chisq + bin2_chisq
    
    return bin1_chisq, bin2_chisq, total_chisq

def stxs_fit(root_name, coffea_name, operator_list, verbose = False):
    start0_MG_sigma = 3.594
    if len(operator_list) > 2:
        raise ValueError(f"Currently only works for 2 or less operators. Your operator list is length: {len(operator_list)}.")
    
    bin1_quad_dependence, bin2_quad_dependence = stxs_reweight_function(coffea_name, operator_list, verbose = verbose)

    h = uproot.open(root_name)["Runs"]
    sumw = h["genEventSumw"].array(library="np").sum()

    wc_space_dict = {}
    
    if len(operator_list) == 2:
        op1 = operator_list[0]
        op2 = operator_list[1]

        for i1 in np.linspace(-10,10,50):
            for i2 in np.linspace(-10,10,50):
                wcpoint = (i1,i2)
                bin1_xsec = start0_MG_sigma * 1000 * quad_2D(wcpoint, *bin1_quad_dependence) / sumw
                bin2_xsec = start0_MG_sigma * 1000 * quad_2D(wcpoint, *bin2_quad_dependence) / sumw
                total_xsec = bin1_xsec + bin2_xsec
    
                bin1_chisq, bin2_chisq, total_chisq = chisq_stxs(bin1_xsec, bin2_xsec)
                
                wc_space_dict[wcpoint] = {'bin1': [bin1_xsec, bin1_chisq],
                                          'bin2': [bin2_xsec, bin2_chisq],
                                          'total': [total_xsec, total_chisq]
                                         }

    elif len(operator_list) == 1:
        op1 = operator_list[0]

        for i1 in np.linspace(-10,10,50):
            wcpoint = (i1)
            bin1_xsec = start0_MG_sigma * 1000 * quad_1D(wcpoint, *bin1_quad_dependence) / sumw
            bin2_xsec = start0_MG_sigma * 1000 * quad_1D(wcpoint, *bin2_quad_dependence) / sumw
            total_xsec = bin1_xsec + bin2_xsec

            bin1_chisq, bin2_chisq, total_chisq = chisq_stxs(bin1_xsec, bin2_xsec)
            
            wc_space_dict[wcpoint] = {'bin1': [bin1_xsec, bin1_chisq],
                                      'bin2': [bin2_xsec, bin2_chisq],
                                      'total': [total_xsec, total_chisq]
                                     }

    best_point = None
    best_chisq = float('inf')
    
    for wcpoint, values in wc_space_dict.items():
        chisq = values['total'][1] 
        if chisq < best_chisq:
            best_chisq = chisq
            best_point = wcpoint
    
    print("Best wc point:", best_point)
    print("Minimum chi2:", best_chisq)
    
    return wc_space_dict


def one_sigma_1d(wc_space_dict, delta_chi2=1.0):
    points = sorted(wc_space_dict.keys())
    chi2_vals = np.array([wc_space_dict[p]['total'][1] for p in points])

    chi2_min = np.min(chi2_vals)
    best_point = points[np.argmin(chi2_vals)]

    allowed_points = [p for p in points if wc_space_dict[p]['total'][1] <= chi2_min + delta_chi2]

    if not allowed_points:
        low = best_point
        high = best_point

    if allowed_points:
        low = min(allowed_points)
        high = max(allowed_points)

    return best_point, low, high

def operator_sigma_1d(root_path, coffea_name, operator_list):
    results = {}

    for op in operator_list:
        wc_space_dict = stxs_fit(root_path, coffea_name, [op])
        best, low, high = one_sigma_1d(wc_space_dict)

        results[op] = [best, low, high]

    return results