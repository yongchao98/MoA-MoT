import numpy as np

# Constants
N_A = 6.0221409e23  # mol-1
R = 8.3144598       # J/(K.mol)
g0 = 9.80665        # m/s2

M_air = 28.965338 * 1e-3  # kg/mol
M_air_dry = 28.9644 * 1e-3  # kg/mol
M_h2o = 18.016 * 1e-3  # kg/mol

P0 = 101325.  # Pa
T0 = 288.15  # K
L = 0.0065  # K/m

def Pressure_isothermal(altitude):
    h = altitude
    P = P0 * np.exp(-((g0 * M_air_dry) / (R * T0)) * h)
    return P

def Pressure_adiabatic(altitude):
    h = altitude
    P = P0 * np.exp(g0 * M_air_dry / R / L * np.log(1 - L * h / T0))
    return P

def XDepth_isothermal(altitude, costh=1):
    h = altitude
    XD = Pressure_isothermal(h) / g0 / costh
    return XD

def XDepth_adiabatic(altitude, costh=1):
    h = altitude
    XD = Pressure_adiabatic(h) / g0 / costh
    return XD

def RayOptDepth_adiabatic(wavelength, altitude=0, costh=1):
    h = altitude
    A = XDepth_adiabatic(h, costh) / (3102. * 1e-3 / (1e-4))
    B = np.exp(-4. * np.log(wavelength / 400.))
    C = 1 - 0.0722 * np.exp(-2 * np.log(wavelength / 400.))
    OD = A * B / C
    return OD

def RayOptDepth_isothermal(wavelength, altitude=0, costh=1):
    h = altitude
    A = XDepth_isothermal(h, costh) / (31020.)
    B = np.exp(-4. * np.log(wavelength / 400.))
    C = 1 - 0.0722 * np.exp(-2 * np.log(wavelength / 400.))
    OD = A * B / C
    return OD

def main_solution(altitude, wavelength, costh=1):
    # Convert inputs to appropriate types if necessary
    altitude = float(altitude)
    wavelength = float(wavelength)
    costh = float(costh)
    
    # Calculate Rayleigh optical depth using isothermal model
    opt_depth_isothermal = RayOptDepth_isothermal(wavelength, altitude, costh)
    
    # Calculate Rayleigh optical depth using adiabatic model
    opt_depth_adiabatic = RayOptDepth_adiabatic(wavelength, altitude, costh)
    
    # Return the results as a dictionary
    return {
        "opt_depth_isothermal": float(opt_depth_isothermal),
        "opt_depth_adiabatic": float(opt_depth_adiabatic)
    }

# Target optical depths
target_isothermal = 0.12503162367592735
target_adiabatic = 0.12378936506099655

# Search for feasible inputs
feasible_altitude = None
feasible_wavelength = None
costh = 1

# Iterate over a wider range of altitudes and wavelengths
for altitude in np.arange(0, 30000, 10):  # Altitude in meters
    for wavelength in np.arange(300, 800, 0.1):  # Wavelength in nanometers
        result = main_solution(altitude, wavelength, costh)
        if (np.isclose(result['opt_depth_isothermal'], target_isothermal, atol=1e-4) and
            np.isclose(result['opt_depth_adiabatic'], target_adiabatic, atol=1e-4)):
            feasible_altitude = altitude
            feasible_wavelength = wavelength
            break
    if feasible_altitude is not None:
        break

print({"altitude": feasible_altitude, "wavelength": feasible_wavelength, "costh": costh})