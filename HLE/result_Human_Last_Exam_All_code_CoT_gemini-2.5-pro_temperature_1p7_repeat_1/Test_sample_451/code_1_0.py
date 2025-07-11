import scipy.constants as const

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for rotational transitions in H-Br.
    """
    # --- Step 1: Define Constants ---
    # Given in the problem
    r0_pm = 141.4  # bond length in picometers
    k_Nm = 400.0   # force constant in N/m
    mH_amu = 1.008 # mass of Hydrogen in amu
    mBr_amu = 79.904 # mass of Bromine in amu

    # Physical constants from scipy.constants for precision
    amu_to_kg = const.physical_constants["atomic mass unit-kilogram relationship"][0]
    h = const.h        # Planck constant in J.s
    c = const.c        # Speed of light in m/s
    e_charge = const.e # Elementary charge in C

    # Unit conversions
    pm_to_m = 1e-12
    eV_to_qeV = 1e30 # Factor to convert eV to qeV (1 qeV = 1e-30 eV)

    # --- Step 2: Calculate Molecular Properties ---
    # Convert inputs to SI units
    r0 = r0_pm * pm_to_m
    mH_kg = mH_amu * amu_to_kg
    mBr_kg = mBr_amu * amu_to_kg

    # Reduced mass (mu) in kg
    mu = (mH_kg * mBr_kg) / (mH_kg + mBr_kg)

    # Moment of inertia (I) in kg.m^2
    I = mu * r0**2

    # --- Step 3: Calculate Spectroscopic Constants (in m^-1)---
    # Rotational constant B
    B_wavenumber = h / (8 * const.pi**2 * I * c)

    # Vibrational frequency nu_e
    omega_e = (k_Nm / mu)**0.5
    nu_e_wavenumber = omega_e / (2 * const.pi * c)

    # Centrifugal distortion constant D
    D_wavenumber = 4 * B_wavenumber**3 / nu_e_wavenumber**2

    # --- Step 4: Calculate Energy Shift for Transitions ---
    # The energy shift formula: ΔE_shift = -4 * h * c * D * (J+1)^3
    
    print("This script calculates the energy shifts due to centrifugal distortion for H-Br.")
    print("Formula for the shift in transition energy: ΔE = -4 * h * c * D * (J+1)^3\n")
    print(f"Derived Spectroscopic Constants for H-Br:")
    print(f"h (Planck constant) = {h:.6e} J·s")
    print(f"c (speed of light)  = {c:.6e} m/s")
    print(f"D (distortion const)= {D_wavenumber:.6e} m^-1\n")
    
    # Transition 1: J = 0 to J = 1
    J0 = 0
    delta_E_J0_to_J1_joules = -4 * h * c * D_wavenumber * (J0 + 1)**3
    delta_E_J0_to_J1_qeV = (delta_E_J0_to_J1_joules / e_charge) * eV_to_qeV

    print("1. Energy shift for transition from J = 0 to J = 1:")
    print(f"   ΔE = -4 * ({h:.6e}) * ({c:.6e}) * ({D_wavenumber:.6e}) * ({J0} + 1)^3")
    print(f"   Final Answer: {delta_E_J0_to_J1_qeV:.4e} qeV\n")

    # Transition 2: J = 1 to J = 2
    J1 = 1
    delta_E_J1_to_J2_joules = -4 * h * c * D_wavenumber * (J1 + 1)**3
    delta_E_J1_to_J2_qeV = (delta_E_J1_to_J2_joules / e_charge) * eV_to_qeV
    
    print("2. Energy shift for transition from J = 1 to J = 2:")
    print(f"   ΔE = -4 * ({h:.6e}) * ({c:.6e}) * ({D_wavenumber:.6e}) * ({J1} + 1)^3")
    print(f"   Final Answer: {delta_E_J1_to_J2_qeV:.4e} qeV")

calculate_energy_shifts()
# The final answer block is below, as requested by the user prompt.
# I'm providing the final computed values directly.
# For J=0->1: -1.7720e+23 qeV
# For J=1->2: -1.4176e+24 qeV