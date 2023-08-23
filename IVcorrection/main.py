import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from typing import Union
from pvlib.pvsystem import calcparams_desoto, singlediode
from pvlib.ivtools.sdm import fit_desoto
from pvlib.ivtools.sde import fit_sandia_simple
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections.abc import Sized

"""
The current version works for the correction of I-V curves of c-Si PV modules. 
The performance on other PV technologies will be updated later.

"""

irradiance_STC = 1000
temperature_STC = 25

colors = ["lightseagreen", "lawngreen", "gold", "darkorange"]
cmap = LinearSegmentedColormap.from_list("mycmap", colors)

    
def simu_IV_curve(G : Union[list, float], 
                T : Union[list, float],
                alpha_isc_abs : float,
                SDMparams : dict,
                ivcurve_pnts : int =100,
                singlediode_method : str ='lambertw'):
    
    """
    Simulate I-V curves at specified environmental conditions

    Parameters
    ----------
    G : float or list
        Irradiance [W/m2]
    T : float or list
        Module temperature [C]
    alpha_isc_abs : float
        Absolute temperature coefficient of Isc (A/C)
    SDMparams: dict
        Dict including the five single-diode model (SDM) parameters
        The five parameters (keys of the dict) are:
            I_L_ref: Light-generated current at reference conditions [A]
            I_o_ref: Diode saturation current at reference conditions [A]
            R_s: Series resistance [ohm]
            R_sh_ref: Shunt resistance at reference conditions [ohm]
            a_ref: Modified ideality factor at reference conditions.
                    The product of the usual diode ideality factor (n, unitless),
                    number of cells in series (Ns), and cell thermal voltage at
                    specified effective irradiance and cell temperature.
    ivcurve_pnts : int 
        Number of points in the desired IV curve. If None or 0, no IV curves will be produced.
    singlediode_method : str
        Determines the method used to calculate points on the IV curve. 
        The options are 'lambertw', 'newton', or 'brentq'.

    Output
    ------
    ivcurves: dict
        Dict including the generated 'v', 'i' and 'G', 'T'

    """

    G = np.array(G)
    T = np.array(T)

    Eg_ref = 1.121
    dEgdT = -0.0002677
    
    iph, io, rs, rsh, nNsVth = calcparams_desoto(
                                G,
                                T,
                                alpha_sc = alpha_isc_abs,
                                a_ref = SDMparams['a_ref'],
                                I_L_ref = SDMparams['I_L_ref'],
                                I_o_ref = SDMparams['I_o_ref'],
                                R_sh_ref = SDMparams['R_sh_ref'],
                                R_s = SDMparams['R_s'],
                                EgRef=Eg_ref,
                                dEgdT=dEgdT,
                                irrad_ref=irradiance_STC,
                                temp_ref=temperature_STC
                                )

    if singlediode_method in ['lambertw', 'brentq', 'newton']:
        out = singlediode(iph,
                        io,
                        rs,
                        rsh,
                        nNsVth,
                        method=singlediode_method,
                        ivcurve_pnts=ivcurve_pnts,
                        )
    else:
        raise Exception(
            'Method must be "lambertw", "brentq", or "newton"')

    ivcurves = {}
    alli = {}
    allv = {}

    nIV = G.size

    if nIV > 1:
        for i in range(nIV):
            alli[i] = out['i'][i]
            allv[i] = out['v'][i]

    elif nIV == 1:
        alli = out['i']
        allv = out['v']
        

    ivcurves['G'] = G
    ivcurves['T'] = T
    ivcurves['v'] = allv
    ivcurves['i'] = alli
    
    return ivcurves

def get_corrected_IV_P1 (iv_initial : dict, 
                    alpha_isc_abs : float,
                    beta_voc_abs : float,
                    rs : float = 0.35,
                    k : float = 0
                    ):
    
    """
    Get the corrected I-V curves using Procedure 1 of IEC 60891:2021 [1]

    Parameters
    ----------
    iv_initial : dict
        Dict including the 'v', 'i' of the I-V curves to correct, 
        and G, T where the I-V curve is measured under
    alpha_isc_abs : float
        Absolute temperature coefficient of short circuit current [A/C]
    beta_voc_abs : float
        Absolute temperature coefficient of open circuit voltage [V/C]
    rs : float 
        Correction coefficient of the internal series resistance [ohm]
    k : float 
        Correction coefficient of the curve correction factor
    
    Output
    ------
    iv_corrected: dict
        Dict including the 'v', 'i' of the corrected I-V curves, 
        and 'G', 'T' where the initial I-V curve is measured under

    Reference
    ---------
    [1] IEC 60891:2021, Photovoltaic devices - Procedures for temperature 
    and irradiance corrections to measured I-V characteristics
    
    """
    
    iv_corrected = {key: iv_initial[key] for key in ['G', 'T']}
    alli_corr = {}
    allv_corr = {}

    nIV = np.array(iv_initial['G']).size
                    
    for n in range(nIV):
        if nIV > 1:
            i = iv_initial['i'][n]
            v = iv_initial['v'][n]
            G = iv_initial['G'][n]
            T = iv_initial['T'][n]
        else:
            i = iv_initial['i']
            v = iv_initial['v']
            G = iv_initial['G']
            T = iv_initial['T']

        isc = np.max(i)

        i_corr = i+isc*(1000/G-1) + alpha_isc_abs*(25-T)
        v_corr = v - rs*(i_corr-i) - k*i_corr*(25-T) + beta_voc_abs*(25-T)

        alli_corr[n] = i_corr
        allv_corr[n] = v_corr

    iv_corrected['v'] = allv_corr
    iv_corrected['i'] = alli_corr

    return iv_corrected

def get_corrected_IV_P2 (iv_initial : dict, 
                    alpha_isc_rel : float,
                    beta_voc_rel : float,
                    voc_ref : float,
                    rs : float = 0.35,
                    k : float = 0,
                    B1 : float = 0,
                    B2 : float = 0):
    
    """
    Get the corrected I-V curves using Procedure 2 of IEC 60891:2021 [1]

    Parameters
    ----------
    iv_initial : dict
        Dict including the 'v', 'i' of the I-V curves to correct, 
        and G, T where the I-V curve is measured under
    alpha_isc_rel : float
        Relative temperature coefficient of short circuit current [A/%]
    beta_voc_rel : float
        Relative temperature coefficient of open circuit voltage [V/%]
    voc_ref : float
        Open-circuit voltage at STC [V]
    rs : float 
        Correction coefficient of the internal series resistance [ohm]
    k : float 
        Correction coefficient of the curve correction factor
    B1 : float 
        Correction coefficient of the irradiance correction factor 1 
    B2 : float 
        Correction coefficient of the irradiance correction factor 2
    
    Output
    ------
    iv_corrected: dict
        Dict including the 'v', 'i' of the corrected I-V curves, 
        and 'G', 'T' where the initial I-V curve is measured under

    Reference
    ---------
    [1] IEC 60891:2021, Photovoltaic devices - Procedures for temperature 
    and irradiance corrections to measured I-V characteristics
    
    """
    
    
    iv_corrected = {key: iv_initial[key] for key in ['G', 'T']}
    alli_corr = {}
    allv_corr = {}

    nIV = np.array(iv_initial['G']).size
                    
    for n in range(nIV):

        if nIV > 1:
            i = iv_initial['i'][n]
            v = iv_initial['v'][n]
            G = iv_initial['G'][n]
            T = iv_initial['T'][n]
        else:
            i = iv_initial['i']
            v = iv_initial['v']
            G = iv_initial['G']
            T = iv_initial['T']

        i_corr = i*1000/G/(1 + alpha_isc_rel*(T-25))
        
        rs1 = rs + k*(T-25)
        fG1 = B2*(np.log(1000/G)**2) + B1*np.log(1000/G) + 1

        v_corr = v - rs1*(i_corr-i) - k*i_corr*(25-T) + \
            voc_ref*(beta_voc_rel*(25-T)*fG1 + 1 - 1/fG1)
        
        alli_corr[n] = i_corr
        allv_corr[n] = v_corr

    iv_corrected['v'] = allv_corr
    iv_corrected['i'] = alli_corr

    return iv_corrected

def get_corrected_IV_P4 (iv_initial : dict,
                      alpha_isc_abs : float,
                      N_cells : int,
                      rs : float = 0.35,
                      epsilon : float = 1.232):

    """
    Get the corrected I-V curves using Procedure 4 of IEC 60891:2021 [1]

    Parameters
    ----------
    iv_initial : dict
        Dict including the 'v', 'i' of the I-V curves to correct, 
        and G, T where the I-V curve is measured under
    alpha_isc_abs : float
        Absolute temperature coefficient of short circuit current [A/C]
    N_cells : int
        Number of cells of PV module
    rs : float 
        Correction coefficient of the internal series resistance [ohm]
    epsilon : float 
        Bandgap of the PV module, default for cs-Si module: 1.232


    Output
    ------
    iv_corrected: dict
        Dict including the 'v', 'i' of the corrected I-V curves, 
        and 'G', 'T' where the initial I-V curve is measured under
    
    Reference
    ---------
    [1] IEC 60891:2021, Photovoltaic devices - Procedures for temperature 
    and irradiance corrections to measured I-V characteristics

    """


    iv_corrected = {key: iv_initial[key] for key in ['G', 'T']}
    alli_corr = {}
    allv_corr = {}

    nIV = np.array(iv_initial['G']).size
                    
    for n in range(nIV):

        if nIV > 1:
            i = iv_initial['i'][n]
            v = iv_initial['v'][n]
            G = iv_initial['G'][n]
            T = iv_initial['T'][n]
        else:
            i = iv_initial['i']
            v = iv_initial['v']
            G = iv_initial['G']
            T = iv_initial['T']

        isc = np.max(i)

        try:
            rs = get_rs_sandia(v,i)

        except:
            alli_corr[n] = i_corr
            allv_corr[n] = v_corr
            
        else:
            itemp = i + isc*(1000/G-1) 
            vtemp = v - rs*(itemp-i)
            i_corr = itemp + alpha_isc_abs*(25-T)
            v_corr = vtemp + (25-T)/(T+273.15)*(vtemp-epsilon*N_cells)

            alli_corr[n] = i_corr
            allv_corr[n] = v_corr

    iv_corrected['v'] = allv_corr
    iv_corrected['i'] = alli_corr

    return iv_corrected

def get_corrected_IV_Pdyna (iv_initial : dict[float],  
                        alpha_isc_rel : float,
                        beta_voc_rel : float,
                        voc_ref : float,
                        rs : float = 0.35,
                        k : float = 0,
                        B1 : float = 0,
                        B2 : float = 0):
    
    """
    Get the corrected I-V curves using Procedure dynamic (Pdyna)

    Parameters
    ----------
    iv_initial : dict
        Dict including the 'v', 'i' of the I-V curves to correct, 
        and 'G', 'T' where the I-V curve is measured under
    alpha_isc_rel : float
        Relative temperature coefficient of short circuit current [A/%]
    beta_voc_rel : float
        Relative temperature coefficient of open circuit voltage [V/%]
    voc_ref : float
        Open-circuit voltage at STC [V]
    rs : float 
        Correction coefficient of the internal series resistance [ohm]
        It is calculatd using the sample_simple funcition from the raw I-V curve [1]
    k : float 
        Correction coefficient of the curve correction factor
    B1 : float 
        Correction coefficient of the irradiance correction factor 1 
    B2 : float 
        Correction coefficient of the irradiance correction factor 2
    
    Output
    ------
    iv_corrected: dict
        Dict including the 'v', 'i' of the corrected I-V curves, 
        and 'G', 'T' where the initial I-V curve is measured under

    Reference
    ---------
    [1] C. B. Jones and C. W. Hansen, “Single Diode Parameter Extraction from In-Field 
        Photovoltaic I-V Curves on a Single Board Computer,” Conference Record of the IEEE 
        Photovoltaic Specialists Conference, pp. 382–387, Jun. 2019, doi: 10.1109/PVSC40753.2019.8981330.
    
    """    

    
    iv_corrected = {key: iv_initial[key] for key in ['G', 'T']}

    alli_corr = {}
    allv_corr = {}

    nIV = np.array(iv_initial['G']).size
                    
    for n in range(nIV):
        
        if nIV > 1:
            i = iv_initial['i'][n]
            v = iv_initial['v'][n]
            G = iv_initial['G'][n]
            T = iv_initial['T'][n]
        else:
            i = iv_initial['i']
            v = iv_initial['v']
            G = iv_initial['G']
            T = iv_initial['T']

        try:
            rs = get_rs_sandia(v[i>0],i[i>0]) # dynamic rs
            rsh = get_rsh_sandia(v[i>0],i[i>0]) # dynamic rsh

            if rsh > 100: # neglect rsh when it is large
                rsh = 1e5

        except:
            alli_corr[n] = np.nan
            allv_corr[n] = np.nan
            
        else:
            itemp = (i+v/rsh)*1000/G/(1 + alpha_isc_rel*(T-25))
            
            rs1 = rs + k*(T-25)
            fG1 = B2*(np.log(1000/G)**2) + B1*np.log(1000/G) + 1

            v_corr = v - rs1*(itemp-i-v/rsh) - k*itemp*(25-T) + \
                voc_ref*(beta_voc_rel*(25-T)*fG1 + 1 - 1/fG1)
            
            i_corr = itemp -v_corr/rsh

            alli_corr[n] = i_corr
            allv_corr[n] = v_corr

    iv_corrected['v'] = allv_corr
    iv_corrected['i'] = alli_corr

    return iv_corrected

def get_rs_coef(SDMparams : dict,
              alpha_isc_abs : float,
              beta_voc_abs : float = 0,
              alpha_isc_rel : float = 0,
              beta_voc_rel : float = 0,
              voc_ref : float = 0,
              B1 : float = 0,
              B2 : float = 0,
              correction_method : str = None):

    """
    Estimate the correction coefficient rs for Procedure 1 and 2 based on IEC 60891:2021 standard [1]

    Parameters
    ----------
    
    SDMparams: dict
        Dict including the five single-diode model (SDM) parameters
        The five parameters (keys of the dict) are:
            I_L_ref: Light-generated current at reference conditions [A]
            I_o_ref: Diode saturation current at reference conditions [A]
            R_s: Series resistance [ohm]
            R_sh_ref: Shunt resistance at reference conditions [ohm]
            a_ref: Modified ideality factor at reference conditions.
                    The product of the usual diode ideality factor (n, unitless),
                    number of cells in series (Ns), and cell thermal voltage at
                    specified effective irradiance and cell temperature.
    alpha_isc_abs : float
        Absolute temperature coefficient of short-circuit current (A/C)
    beta_voc_abs : float
        Absolute temperature coefficient of open circuit voltage [V/C]
    alpha_isc_rel : float
        Relative temperature coefficient of short circuit current [A/%]
    beta_voc_rel : float
        Relative temperature coefficient of open circuit voltage [V/%]
    voc_ref : float
        Open-circuit voltage at STC [V]
    B1 : float 
        Correction coefficient of the irradiance correction factor 1 for Procedure 2
    B2 : float 
        Correction coefficient of the irradiance correction factor 2 for Procedure 2
    correction_method : str
        Correction method name, choose from 'P1' or 'P2'

    Output
    ------
    rs_best: float
        Optimal correction coefficient rs

    Reference
    ---------
    [1] IEC 60891:2021, Photovoltaic devices - Procedures for temperature 
    and irradiance corrections to measured I-V characteristics

    """

    # Simulate IV curves at varying irradiance with temperature = 25
    G = np.arange(200, 1201, 100)
    T = np.ones(len(G))*25

    ivcurves = simu_IV_curve(G, T, alpha_isc_abs, SDMparams)
    ivcurve_STC = simu_IV_curve([irradiance_STC], [temperature_STC], alpha_isc_abs, SDMparams)

    if correction_method is None:
        raise ValueError('Correction method must be "P1" or "P2"')

    rs_range = np.arange(0.1, 0.7, 0.01)

    rs_results = {}
    rs_results['rs'] = rs_range
    rs_results['pmp_err_mean'] = [np.nan] * len(rs_range)
    rs_results['pmp_err_std'] = [np.nan] * len(rs_range)

    for i in range(len(rs_range)):
        if correction_method == 'P1':
            ivcurves_corr = get_corrected_IV_P1(ivcurves, 
                                        alpha_isc_abs = alpha_isc_abs,
                                        beta_voc_abs = beta_voc_abs,
                                        rs = rs_range[i])
        elif correction_method == 'P2':
            ivcurves_corr = get_corrected_IV_P2(ivcurves, 
                                        alpha_isc_rel = alpha_isc_rel,
                                        beta_voc_rel = beta_voc_rel,
                                        voc_ref= voc_ref,
                                        B1 = B1,
                                        B2 = B2,
                                        rs = rs_range[i])
        
        allpmp_err = np.abs(calc_mpp_error(ivcurves_corr, ivcurve_STC))
        rs_results['pmp_err_mean'][i] = np.mean(allpmp_err)
        rs_results['pmp_err_std'][i] = np.std(allpmp_err)

    rs_best = rs_results['rs'][np.argmin(rs_results['pmp_err_mean'])]
    
    return rs_best, rs_results

def get_k_coef(SDMparams : dict,
              alpha_isc_abs : float,
              beta_voc_abs : float = 0,
              alpha_isc_rel : float = 0,
              beta_voc_rel : float = 0,
              voc_ref : float = 0,
              B1 : float = 0, 
              B2 : float = 0,
              rs : float = 0,
              correction_method : str = None):
    
    """
    Estimate the correction coefficient k for Procedure 1 or 2 based on IEC 60891:2021 standard [1]

    Parameters
    ----------
    SDMparams: dict
        Dict including the five single-diode model (SDM) parameters
        The five parameters (keys of the dict) are:
            I_L_ref: Light-generated current at reference conditions [A]
            I_o_ref: Diode saturation current at reference conditions [A]
            R_s: Series resistance [ohm]
            R_sh_ref: Shunt resistance at reference conditions [ohm]
            a_ref: Modified ideality factor at reference conditions.
                    The product of the usual diode ideality factor (n, unitless),
                    number of cells in series (Ns), and cell thermal voltage at
                    specified effective irradiance and cell temperature.
    alpha_isc_abs : float
        Absolute temperature coefficient of short-circuit current (A/C)
    beta_voc_abs : float
        Absolute temperature coefficient of open circuit voltage [V/C]
    alpha_isc_rel : float
        Relative temperature coefficient of short circuit current [A/%]
    beta_voc_rel : float
        Relative temperature coefficient of open circuit voltage [V/%]
    voc_ref : float
        Open-circuit voltage at STC [V]
    B1 : float 
        Correction coefficient of the irradiance correction factor 1 for Procedure 2
    B2 : float 
        Correction coefficient of the irradiance correction factor 2 for Procedure 2
    rs: float
        Correction coefficient rs
    correction_method : str
        Correction method name, choose from 'P1' or 'P2'

    Output
    ------
    k_best: float
        Optimal correction coefficient k

    Reference
    ---------
    [1] IEC 60891:2021, Photovoltaic devices - Procedures for temperature 
    and irradiance corrections to measured I-V characteristics

    """          

    # Simulate IV curves at varying temperature with irradiance = 1000W/m2
    T = np.arange(10, 71, 5)
    G = np.ones(len(T))*1000

    ivcurves = simu_IV_curve(G, T, alpha_isc_abs, SDMparams)
    ivcurve_STC = simu_IV_curve([irradiance_STC], [temperature_STC], alpha_isc_abs, SDMparams)

    if correction_method is None:
        raise ValueError('Correction method must be "P1" or "P2"')
    
    k_range = np.arange(0, 0.01, 0.0001)

    k_results = {}
    k_results['k'] = k_range
    k_results['pmp_err_mean'] = [np.nan] * len(k_range)
    k_results['pmp_err_std'] = [np.nan] * len(k_range)

    for i in range(len(k_range)):
        if correction_method == 'P1':
            ivcurves_corr = get_corrected_IV_P1(ivcurves, 
                                    alpha_isc_abs = alpha_isc_abs,
                                    beta_voc_abs = beta_voc_abs,
                                    rs = rs, 
                                    k = k_range[i])
            
        elif correction_method == 'P2':
            ivcurves_corr = get_corrected_IV_P2(ivcurves, 
                                        alpha_isc_rel = alpha_isc_rel,
                                        beta_voc_rel = beta_voc_rel,
                                        voc_ref= voc_ref,
                                        B1 = B1,
                                        B2 = B2,
                                        rs = rs,
                                        k = k_range[i])

        allpmp_err = np.abs(calc_mpp_error(ivcurves_corr, ivcurve_STC))
        k_results['pmp_err_mean'][i] = np.mean(allpmp_err)
        k_results['pmp_err_std'][i] = np.std(allpmp_err)

    k_best = k_range[np.argmin(k_results['pmp_err_std'])]

    return k_best, k_results

def get_rs_sandia(v : list, i : list):

    """
    Estimate the correction coefficient rs from the raw I-V curve using sandia_simple [1]

    Parameters
    ----------
    v : list 
        Voltage of I-V curve [V]
    i : list 
        Current of I-V curve [A]

    Output
    ------
    rs: float
        Correction coefficient rs [ohms]

    Reference
    ---------
    [1] C. B. Jones and C. W. Hansen, “Single Diode Parameter Extraction from In-Field 
        Photovoltaic I-V Curves on a Single Board Computer,” Conference Record of the IEEE 
        Photovoltaic Specialists Conference, pp. 382–387, Jun. 2019, doi: 10.1109/PVSC40753.2019.8981330.

    """    

    rs = fit_sandia_simple(v, i)[2]
    return rs

def get_rsh_sandia(v : list, i : list):

    """
    Estimate the correction coefficient rsh from the raw I-V curve using sandia_simple [1]

    Parameters
    ----------
    v : list 
        Voltage of I-V curve [V]
    i : list 
        Current of I-V curve [A]

    Output
    ------
    rsh: float
        Correction coefficient rsh [ohms] 

    Reference
    ---------
    [1] C. B. Jones and C. W. Hansen, “Single Diode Parameter Extraction from In-Field 
        Photovoltaic I-V Curves on a Single Board Computer,” Conference Record of the IEEE 
        Photovoltaic Specialists Conference, pp. 382–387, Jun. 2019, doi: 10.1109/PVSC40753.2019.8981330.

    """ 

    rsh = fit_sandia_simple(v, i)[3]
    return rsh

def get_B1B2_coef(alpha_isc_abs : float,
                    voc_ref : float,
                    SDMparams : dict):

    """
    Estimate the correction coefficient B1 and B2 for Procedure 2 based on IEC 60891:2021 standard [1]

    Parameters
    ----------
    alpha_isc_abs : float
        Absolute temperature coefficient of short-circuit current (A/C)
    voc_ref : float
        Open-circuit voltage at STC [V]
    SDMparams: dict
        Dict including the five single-diode model (SDM) parameters
        The five parameters (keys of the dict) are:
            I_L_ref: Light-generated current at reference conditions [A]
            I_o_ref: Diode saturation current at reference conditions [A]
            R_s: Series resistance [ohm]
            R_sh_ref: Shunt resistance at reference conditions [ohm]
            a_ref: Modified ideality factor at reference conditions.
                    The product of the usual diode ideality factor (n, unitless),
                    number of cells in series (Ns), and cell thermal voltage at
                    specified effective irradiance and cell temperature.

    Output
    ------
    B1 : float 
        Correction coefficient of the irradiance correction factor 1 
    B2 : float 
        Correction coefficient of the irradiance correction factor 2

    Reference
    ---------
    [1] IEC 60891:2021, Photovoltaic devices - Procedures for temperature 
    and irradiance corrections to measured I-V characteristics

    """  

    # Simulate IV curves at varying irradiance with temperature = 25
    G = np.arange(400, 1101, 100)
    T = np.ones(len(G))*25

    ivcurves = simu_IV_curve(G, T, alpha_isc_abs, SDMparams)

    allvoc = [np.max(ivcurves['v'][n]) for n in range(len(G))]

    x = np.log(1000/G)
    y = voc_ref/np.array(allvoc)

    para = np.polyfit(x, y, 2)

    B1 = para[1]
    B2 = para[0]

    return B1, B2

def get_P1_coefs(alpha_isc_abs : float,
                 beta_voc_abs : float,
                 SDMparams : dict):

    """
    Run pipeline to get rs and k for Procedure 1

    Parameters
    ----------
    alpha_isc_abs : float
        Absolute temperature coefficient of short-circuit current (A/C)
    beta_voc_abs : float
        Absolute temperature coefficient of open circuit voltage [V/C]
    SDMparams: dict
        Dict including the five single-diode model (SDM) parameters
        The five parameters (keys of the dict) are:
            I_L_ref: Light-generated current at reference conditions [A]
            I_o_ref: Diode saturation current at reference conditions [A]
            R_s: Series resistance [ohm]
            R_sh_ref: Shunt resistance at reference conditions [ohm]
            a_ref: Modified ideality factor at reference conditions.
                    The product of the usual diode ideality factor (n, unitless),
                    number of cells in series (Ns), and cell thermal voltage at
                    specified effective irradiance and cell temperature.

    Output
    ------
    coeff : dict
        Dict including rs and k

    """  

    coeff = {}
    coeff['rs'], _ = get_rs_coef(SDMparams = SDMparams,
                            alpha_isc_abs = alpha_isc_abs, 
                            beta_voc_abs = beta_voc_abs, 
                            correction_method= 'P1')
    coeff['k'], _ = get_k_coef(SDMparams = SDMparams,
                            alpha_isc_abs = alpha_isc_abs, 
                            beta_voc_abs = beta_voc_abs, 
                            rs = coeff['rs'], 
                            correction_method='P1')
    return coeff

def get_P2_coefs(alpha_isc_abs : float,
                 alpha_isc_rel : float,
                 beta_voc_rel : float,
                 voc_ref : float,
                 SDMparams : dict):

    """
    Run pipeline to get B1, B2, rs and k for Procedure 2

    Parameters
    ----------
    alpha_isc_abs : float
        Absolute temperature coefficient of short-circuit current (A/C)
    alpha_isc_rel : float
        Relative temperature coefficient of short circuit current [A/%]
    beta_voc_rel : float
        Relative temperature coefficient of open circuit voltage [V/%]
    voc_ref : float
        Open-circuit voltage at STC [V]
    SDMparams: dict
        Dict including the five single-diode model (SDM) parameters
        The five parameters (keys of the dict) are:
            I_L_ref: Light-generated current at reference conditions [A]
            I_o_ref: Diode saturation current at reference conditions [A]
            R_s: Series resistance [ohm]
            R_sh_ref: Shunt resistance at reference conditions [ohm]
            a_ref: Modified ideality factor at reference conditions.
                    The product of the usual diode ideality factor (n, unitless),
                    number of cells in series (Ns), and cell thermal voltage at
                    specified effective irradiance and cell temperature.
    
    Output
    ------
    coeff : dict
        Dict including B1, B2, rs and k

    """  

    coeff = {}
    coeff['B1'], coeff['B2'] = get_B1B2_coef(SDMparams = SDMparams,
                                                alpha_isc_abs = alpha_isc_abs,
                                                voc_ref = voc_ref)
    coeff['rs'], _ = get_rs_coef(SDMparams = SDMparams,
                         alpha_isc_abs = alpha_isc_abs,
                         alpha_isc_rel = alpha_isc_rel, 
                         beta_voc_rel = beta_voc_rel,
                         voc_ref = voc_ref,
                         B1 = coeff['B1'], 
                         B2 = coeff['B2'],
                         correction_method = 'P2')
    coeff['k'], _ = get_k_coef(SDMparams = SDMparams,
                            alpha_isc_abs = alpha_isc_abs,
                            alpha_isc_rel = alpha_isc_rel, 
                            beta_voc_rel = beta_voc_rel,
                            voc_ref = voc_ref,
                            B1 = coeff['B1'], 
                            B2 = coeff['B2'], 
                            rs = coeff['rs'], 
                            correction_method='P2')
    return coeff

def calc_mpp_error(ivcurves : dict,
                   ivcurve_ref : dict,
                   method : str = 'optimized'):
    """
    Calculate the maximum power point error 

    Parameters
    ----------
    ivcurves : dict
        Dict including the 'v', 'i' of the corrected I-V curves and 'G', 'T'
    ivcurve_ref : dict
        Dict including the 'v', 'i' of the reference I-V curve and 'G', 'T'
    method: str (default: 'optimized')
        'standard': calculate the pmp error
        'optimized': calculate the RMSE of vmp and imp error, 
                    which considers both vmp and imp and can make the 
                    maximum power point of corrected curves closer to 
                    the reference one

    Output
    ------
    allpmp_err : list 
        Ratio of the pmp error compared to the reference pmp [%]

    """

    allmpp_err = []
    i_ref = ivcurve_ref['i'][0]
    v_ref = ivcurve_ref['v'][0]

    for n in range(len(ivcurves['G'])):
        i = ivcurves['i'][n]
        v = ivcurves['v'][n]

        if method == 'standard':
            pmp = np.max(v*i)
            pmp_ref = np.max(ivcurve_ref['v'][0]*ivcurve_ref['i'][0])
            mpp_err = np.abs((pmp-pmp_ref)/pmp_ref*100)

        if method == 'optimized':
            idpmp = np.argmax(v*i)
            idpmp_ref = np.argmax(v_ref*i_ref)
            vmp_err = (v[idpmp]-v_ref[idpmp_ref])/v_ref[idpmp_ref]
            imp_err = (i[idpmp]-i_ref[idpmp_ref])/i_ref[idpmp_ref]
            mpp_err = np.sqrt((vmp_err**2+imp_err**2)/2)

        allmpp_err.append(mpp_err)

    return allmpp_err

def calc_correction_error(ivcurves : dict,
                          iv_ref: dict):
    
    """
    Calculate the area error between the corrected and the reference I-V curve 
    Area error is a ratio of the area size of the difference to the area size of the reference curve [1]

    Parameters
    ----------
    ivcurves : dict
        Dict including the 'v', 'i' of the corrected I-V curves and 'G', 'T'
    ivcurve_ref : dict
        Dict including the 'v', 'i' of the reference I-V curve and 'G', 'T'

    Output
    ------
    area_err : list 
        List of area errors of the corrected I-V curves [%]

    Reference
    ---------
    [1] J. C. H. Phang and D. S. H. Chan, “A review of curve fitting error criteria 
        for solar cell I-V characteristics,” Solar Cells, vol. 18, no. 1, pp. 1–12, 
        Jul. 1986, doi: 10.1016/0379-6787(86)90002-5.

    """

    area_err = []
    nIV = np.array(ivcurves['G']).size

    for n in range(nIV):
        i = ivcurves['i'][n]
        v = ivcurves['v'][n]
        i_ref = iv_ref['i'][0]
        v_ref = iv_ref['v'][0]

        vmin = max(0, np.min(v))

        if np.min(i) > 0:
            vreftemp = np.interp(np.min(i), np.flip(i_ref), np.flip(v_ref))
            vmax = max(np.max(v), vreftemp)
        else:
            voc = np.interp(0, i, v)
            vmax = max(voc, max(v_ref))

        vi = np.arange(vmin, vmax+1e-3, (vmax-vmin)/1000)
        ii = np.interp(vi, v, i)
        irefi = np.interp(vi, v_ref, i_ref)

        icut = max(0, min(i))
        irefi[irefi<icut] = icut
        ii[ii<icut] = icut

        area_err.append(np.sum(np.abs(ii-irefi))/np.sum(irefi-icut)*100)

    return area_err

def plot_iv_corrected (ivcurves : dict, 
                       iv_ref : dict = None, 
                       lw : float = 2, 
                       plot_ref : bool = True, 
                       xmax : float = None, 
                       ymax : float = None, 
                       errorshow : bool = True):
    
    """
    Plot the corrected I-V curves

    Parameters
    ----------
    ivcurves : dict
        Dict including the 'v', 'i' of the corrected I-V curves and 'G', 'T'
    iv_ref : dict
        Dict including the 'v', 'i' of the reference I-V curve and 'G', 'T'
    lw : float
        Line width of the corrected I-V curves
    plot_ref : bool
        Plot the reference curve if True
    xmax : float
        Max of the x-axis of figure
    ymax : float
        Max of the y-axis of figure
    errorshow : bool
        Calculate and display the correction error if True

    """


    plt.rcParams.update({'font.size': 15})
    _, ax = plt.subplots(figsize = [6,4])
    # ax.grid()

    allG = np.array(ivcurves['G'])
    nIV = allG.size
    
    for k in range(nIV):

        color = cmap(k / nIV)  # Get color from colormap

        i = ivcurves['i'][k]
        v = ivcurves['v'][k]

        if k == 0:
            ax.plot(v,i, linewidth = lw, label = 'Corrected', color = color)
        else:
            ax.plot(v,i, linewidth = lw, color = color)

    if xmax:
        ax.set_xlim([0, xmax])
    if ymax:
        ax.set_ylim([0, ymax])

    # if iv_ref is not None:
    #     v_ref = iv_ref['v'][0]
    #     i_ref = iv_ref['i'][0]
    #     if plot_ref:
    #         ax.plot(v_ref, i_ref, '--', color = 'black', linewidth = 1.5, label = 'Ref')
    #     if errorshow:
    #         area_err = calc_correction_error(ivcurves, iv_ref)
    #         ax.text(1, 3, 'Area error: {:.2f}%'.format(np.nanmean(area_err)))

    # Show colorbar when more than 1 curves are corrected
    if nIV > 1:
        # Calculate colorbar tick values based on specific values associated with lines
        colorbar_ticks = np.arange(0, 1.01, 1/4)
        G_bar = np.arange(np.min(allG), np.max(allG)+1, (np.max(allG)-np.min(allG))/4)
        colorbar_tick_labels = [f'{i:.0f}' for i in G_bar]


        cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap), 
                            ax=ax, 
                            ticks=colorbar_ticks,
                            label='Irradiance ($\mathregular{W/m^2}$)')    
        cbar.ax.set_yticklabels(colorbar_tick_labels)

    ax.set_xlim(left = 0)
    ax.set_ylim(bottom = 0)
    ax.set_xlabel('Voltage (V)')
    ax.set_ylabel('Current (A)')
    ax.set_title('I-V curves', fontweight = 'bold')
    # ax.legend(loc=3)  
        
def plot_iv_corrected_quad (iv1 : dict, 
                            iv2 : dict, 
                            iv4 : dict, 
                            ivdyna : dict, 
                            iv_ref : dict = None, 
                            lw : float = 2, 
                            plot_ref : bool = True, 
                            errorshow : bool = True,
                            xmax : float = None,  
                            ymax : float = None):
    

    """
    Plot the corrected I-V curve using 4 correction methods in individual figures

    Parameters
    ----------
    iv1 : dict
        Dict including the 'v', 'i' of the corrected I-V curves using Procedure 1 and 'G', 'T'
    iv2 : dict
        Dict including the 'v', 'i' of the corrected I-V curves using Procedure 2 and 'G', 'T'
    iv4 : dict
        Dict including the 'v', 'i' of the corrected I-V curves using Procedure 4 and 'G', 'T'
    ivdyna : dict
        Dict including the 'v', 'i' of the corrected I-V curves using Pdynamic and 'G', 'T'
    iv_ref : dict
        Dict including the 'v', 'i' of the reference I-V curve and 'G', 'T'
    lw : float
        Line width of the corrected I-V curves
    plot_ref : bool
        Plot the reference curve if True
    xmax : float
        Max of the x-axis of figure
    ymax : float
        Max of the y-axis of figure
    errorshow : bool
        Calculate and display the correction error if True

    """
    
    plt.rcParams.update({'font.size': 11})
    fig, allax = plt.subplots(1, 4, figsize=[16,2.5])
    plt.subplots_adjust(wspace=0.4, hspace=0.1)

    methods = ['P1', 'P2', 'P4', 'Pdyna']

    for n in range(4):
        ax = allax[n]
        methodname = methods[n]

        if n == 0:
            ivcurves = iv1
        elif n == 1:
            ivcurves = iv2
        elif n == 2:
            ivcurves = iv4
        else:
            ivcurves = ivdyna
        
        allG = np.array(ivcurves['G'])
        nIV = allG.size

        for k in range(nIV):
        
            i = ivcurves['i'][k]
            v = ivcurves['v'][k]

            color = cmap(k / nIV)  # Get color from colormap

            if k == 0:
                ax.plot(v,i, linewidth = lw, label = 'Corrected', color = color)
            else:
                ax.plot(v,i, linewidth = lw, color = color)
            if xmax:
                ax.set_xlim([0, xmax])
            if ymax:
                ax.set_ylim([0, ymax])

        if iv_ref is not None:
            if plot_ref:
                ax.plot(iv_ref['v'][0], iv_ref['i'][0], '--', color = 'black', linewidth = 1.5, label = 'Ref')
            if errorshow:
                area_err = calc_correction_error(ivcurves, iv_ref)
                ax.text(1, 3.5, 'Area error: {:.2f}%'.format(np.nanmean(area_err)))

        # Show colorbar when more than 1 curves are corrected
        if (nIV > 1) & (n==3) :
            # Calculate colorbar tick values based on specific values associated with lines
            colorbar_ticks = np.arange(0, 1.01, 1/4)
            G_bar = np.arange(np.min(allG), np.max(allG)+1, (np.max(allG)-np.min(allG))/4)
            colorbar_tick_labels = [f'{i:.0f}' for i in G_bar]

            # Create divider for existing axes instance
            divider = make_axes_locatable(allax[-1])

            # Append a new axes on the right side of the original axes
            cax = divider.append_axes("right", size="5%", pad=0.1)

            cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap), 
                                cax=cax, 
                                ticks=colorbar_ticks,
                                label='Irradiance ($\mathregular{W/m^2}$)')    
            cbar.ax.set_yticklabels(colorbar_tick_labels)

        ax.set_xlabel('Voltage (V)')
        ax.set_ylabel('Current (A)')
        ax.set_title('{}'.format( methodname), fontweight = 'bold')
        ax.set_xlim(left = 0)
        ax.grid()
        ax.legend(loc=3)  
        plt.gcf().set_dpi(120)

def plot_iv_raw(ivcurves : dict, 
                lw : float = 2, 
                title : str = None):
    
    """
    Plot the I-V curves to correct

    Parameters
    ----------
    ivcurves : dict
        Dict including the 'v', 'i' of the curves and 'G', 'T'
    lw : float
        Line width of the I-V curves
    title : str
        Title of the figure

    """

    plt.rcParams.update({'font.size': 15})
    _, ax = plt.subplots(figsize = [6,4])
    # ax.grid()

    allG = np.array(ivcurves['G'])
    nIV = allG.size

    for k in range(nIV):

        color = cmap(k / nIV)  # Get color from colormap

        if nIV == 1:
            i = ivcurves['i']
            v = ivcurves['v']
            G = ivcurves['G']
        else: 
            i = ivcurves['i'][k]
            v = ivcurves['v'][k]
            G = ivcurves['G'][k]

        ax.plot(v, i, color = color, linewidth = lw)
        ax.set_xlim([0, 45])
        ax.set_xlabel('Voltage (V)')
        ax.set_ylabel('Current (A)')
        if title is None:
            ax.set_title('I-V curves', fontweight = 'bold') 
        else:
            ax.set_title(title, fontweight = 'bold')

    if nIV > 1:
        # Calculate colorbar tick values based on specific values associated with lines
        colorbar_ticks = np.arange(0, 1.01, 1/4)
        G_bar = np.arange(np.min(allG), np.max(allG)+1, (np.max(allG)-np.min(allG))/4)
        colorbar_tick_labels = [f'{i:.0f}' for i in G_bar]


        cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap), 
                            ax=ax, 
                            ticks=colorbar_ticks,
                            label='Irradiance ($\mathregular{W/m^2}$)')    
        cbar.ax.set_yticklabels(colorbar_tick_labels)
    