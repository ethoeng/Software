import numpy as np

"""
Required steering angle required to clear-off the aperture,
in units of [mrad],
positive value: steering upwards
negative value: steering downwards
"""
theta_spec_mrad = {'YCB6': 70,
                      'YCB6B': -130-70, #-130 from horizontal axis, -70 to cancel previous divergence
                      'HV6':+130 #steer downward beam to horizontal
                      }

steerer_names = ['YCB6','YCB6B','HV6']

"""
Curent-Field Relation for HH6 from bNMR 1f field calibration measurements in August 2021
(from Derek Fujimoto's bfit HH6 calculator)
"""
def current2field(current):    #in Tesla
    return (current*3.561129+1.574797)*1e-4

def field2current(field):      #in Amps
    return (field-1.574797)/3.561129


# Radial distance of steerer center from sample-2, in [mm]
radialdist_from_sample2 = {'YCB6':497.88,
                           'YCB6B':np.sqrt((497.88-291.88)**2 + 20.41**2),
                           'HV6':(497.88-447.39)}

z_coord_mm = {'YCB6': 497.88, 
              'YCB6B': (497.88-291.88), 
              'HV6': (497.88-447.39),
              'sample2':0}

x_coord_mm = {'YCB6': 0, 
              'YCB6B': 20.41, 
              'HV6': 0,
              'sample2':0}

# x_coord_mm['YCB6B_max'] =
# x_coord_mm['YCB6B_max'] =


"""
Stray magnetic field steering effect -> theta_mag
"""


# Effective length: e.g. int(B,dL)['HV6-YCB6B'] = L_eff*B_at_sample2 
# Positive -> theta_mag upwards, negative -> theta_mag downwards, same convention as theta_spec
def L_eff(steerer_name):
    integration_range = {'YCB6B-end': 7.69675986e-4, 'HV6-YCB6B':-1.39538255e-1, 'sample2-HV6':-5.0603516e-2} 
    if steerer_name=='YCB6': 
        L_eff = integration_range['YCB6B-end']
    elif steerer_name=='YCB6B': 
        L_eff = integration_range['HV6-YCB6B']
    elif steerer_name=='HV6': 
        L_eff = integration_range['sample2-HV6']    
    else:
        raise ValueError(" name not found, available options (with single quote ''): 'YCB6', 'YCB6B','HV6'")
    
    return L_eff

# mag_rigidity = B*rho = sqrt(2*V_0*rest_mass_eV)/c
def mag_rigidity_Teslameter(beam_energy_eV=20e3, isotope_mass_in_nucleon=8):
    return np.sqrt(2*beam_energy_eV*isotope_mass_in_nucleon*931e6)/3e8

# theta_mag = int(B dz)/mag_rigidity
def angle_mag_mrad(HH6_current_A, beam_energy_eV=20e3, isotope_mass_in_nucleon=8):
    """
    Returns dict of divergence w/ keys: 'YCB6','YCB6B','HV6' due to HH6 stray magnetic fields.
    Units in [mrad].
    """
    theta_mag_mrad = dict.fromkeys(steerer_names)
    mag_rigidity_Tm = mag_rigidity_Teslameter(beam_energy_eV, isotope_mass_in_nucleon)

    for name in steerer_names:
        integ_B_dL = current2field(HH6_current_A)*L_eff(name)

        theta_mag_mrad[name] = (integ_B_dL/mag_rigidity_Tm)*1e3

    return theta_mag_mrad 

# steerer dimensions in inches
s = {'YCB6':2,
     'YCB6B':0.875,
     'HV6':1.75} #HV6 gap is estimated as 0.5*(largest_gap_val + smallest_gap_val)

L = {'YCB6':2,
     'YCB6B':1.5,
     'HV6':1}

# V_steer = theta_steerer*(V_0*2s)/L
def get_steerer_voltage_Volt(HH6_current_A, beam_energy_eV=20e3, isotope_mass_in_nucleon=7):
    """
    Returns dict of steerers' gap voltages with dict.keys():'YCB6','YCB6B','HV6'.
    Units of [Volt].
    """
    steerer_voltage = dict.fromkeys(steerer_names)
    steerer_angle = dict.fromkeys(steerer_names)

    theta_mag_mrad = angle_mag_mrad(HH6_current_A, beam_energy_eV, isotope_mass_in_nucleon)

    for name in steerer_names:
        steerer_angle[name] = theta_spec_mrad[name]-theta_mag_mrad[name] 
        steerer_voltage[name] = -(steerer_angle[name]*1e-3)*beam_energy_eV*2*s[name]/L[name]  # positive gap voltage = downward steering (neg angle)
    return steerer_angle, steerer_voltage

def get_magfield_lowerlimit_Tesla(V_gapmax_dict_Volt, beam_energy_eV=20e3, isotope_mass_in_nucleon=8):
    # theta_steerer_mrad = dict.fromkeys(steerer_names)
    # theta_mag_mrad = dict.fromkeys(steerer_names)
    B_min_Tesla = dict.fromkeys(steerer_names)

    for name in steerer_names:
        theta_steerer_mrad = V_gapmax_dict_Volt[name]*L[name]/(beam_energy_eV*2*s[name])
        theta_mag_mrad = theta_spec_mrad[name]- theta_steerer_mrad
        B_min_Tesla[name] = theta_mag_mrad*mag_rigidity_Teslameter(beam_energy_eV, isotope_mass_in_nucleon)/L_eff(name)

    return B_min_Tesla

def get_tot_angle(HH6_current_A, steerer_voltage_dict, beam_energy_eV=20e3, isotope_mass_in_nucleon=8):
    angle_total_mrad = dict.fromkeys(steerer_names)

    theta_mag_mrad = angle_mag_mrad(HH6_current_A, beam_energy_eV, isotope_mass_in_nucleon)
    for key in steerer_voltage_dict.keys():
        angle_total_mrad[key] = theta_mag_mrad[key] + steerer_voltage_dict[key]*L[key]/(beam_energy_eV*2*s[key])
    
    return angle_total_mrad

# def scaleYCB6BfromYCB6(HH6_current_A, YCB6_voltage, beam_energy_eV = 20e3, isotope_mass_in_nucleon=8):
#     V_YCB6B_over_YCB6 = 

#######################################
"""
Run as script for diagnostic purpose
"""
if (__name__ == '__main__'):
    import matplotlib.pyplot as plt
    # from matplotlib impor
    # import argparse

    # user input from stdin
    Beam_erg_eV_input = float(input("Enter beam energy value in eV: "))
    Isotope_mass_nucleon_input = float(input("Enter isotope mass in nucleon: "))
    opt = input("Enter option: \n voltage -> find steerer voltages settings \n limit -> find min HH6 fields possible given steerers gap voltages: ")
    
    if opt == 'voltage':
        # as
        HH6_A_input = current2field(float(input("Enter HH6 magnetic field value in Tesla: ")))

        _, steerer_voltages_V = get_steerer_voltage_Volt(HH6_A_input, Beam_erg_eV_input, Isotope_mass_nucleon_input)
        print("Steerer Voltages are:\n")
        for key, val in zip(steerer_voltages_V.keys(), steerer_voltages_V.values()):
            print(key,':\t', val,' [Volt]\n')

    # elif opt == 'plot':
    #     # GUI for interactive plot here with slider
    #     for key in steerer_names.keys():
    #     """Draw steerers as horizontal lines in x,z planes"""
    #         x0, z0 = zip(x_coord_mm[key], z_coord_mm[key])
    #         x1, z1 = zip(x_coord_mm[key]-L[key]/2,z_coord_mm[key]-s[key]/2) #bottom left corner of steerer plate
    #         x2, z2 = zip(x_coord_mm[key]+L[key]/2,z_coord_mm[key]+s[key]/2) #top right corner of steerer plate
    #         plt.plot(x1, z1, x2, z1) #plot bottom steerer plate
    #         plt.plot(x1, z2, x2, z2) #plot top steerer plate
    #         plt.plot(x0,z0,marker='x') #plot center point of steerer
    
    elif opt == 'limit':
        Vgap_max = dict.fromkeys(steerer_names)
        for name in steerer_names:
            Vgap_max[name] = float(input("Enter max gap voltage available for steerer {} [Volt]: ".format(name)))
        
        B_min = get_magfield_lowerlimit_Tesla(Vgap_max, Beam_erg_eV_input,Isotope_mass_nucleon_input)
        B_min_val  = max(B_min.values()) # find the largest values of magnetic fields amongst B_min[keys]
        print("Minimum fields possible with these gap voltages are: {} mT or {} Gauss.".format(B_min_val*1e-3, B_min_val*1e-2))
    
    elif opt == 'plot':
        # GUI for interactive plot here with slider
        for key in steerer_names:
        # """Draw steerers as horizontal lines in x,z planes"""
            x0, z0 = zip(x_coord_mm[key], z_coord_mm[key])
            x1, z1 = zip(x_coord_mm[key]-L[key]/2,z_coord_mm[key]-s[key]/2) #bottom left corner of steerer plate
            x2, z2 = zip(x_coord_mm[key]+L[key]/2,z_coord_mm[key]+s[key]/2) #top right corner of steerer plate
            plt.plot(x1, z1, x2, z1) #plot bottom steerer plate
            plt.plot(x1, z2, x2, z2) #plot top steerer plate
            plt.plot(x0,z0,marker='x') #plot center point of steerer

    else:
        print("Available options are: voltage, limit")
    


