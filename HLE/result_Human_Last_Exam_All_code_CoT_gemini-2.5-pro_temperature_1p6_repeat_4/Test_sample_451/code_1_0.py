import math
import scipy.constants

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br
    rotational transitions.
    """
    # --- Given Parameters ---
    r0_pm = 141.4      # bond length in pm
    k_force = 400.0    # force constant in N/m
    m_H_amu = 1.008    # mass of hydrogen in amu
    m_Br_amu = 79.904  # mass of bromine in amu

    # --- Physical Constants (from scipy.constants for high precision) ---
    h = scipy.constants.h                         # Planck's constant (J*s)
    pi = scipy.constants.pi                       # Pi
    amu_to_kg = scipy.constants.physical_constants['unified atomic mass unit'][0]
    pm_to_m = 1e-12                               # picometers to meters
    J_to_eV = 1 / scipy.constants.e               # Conversion from Joules to eV
    eV_to_qeV = 1e30                              # Conversion from eV to qeV
    J_to_qeV = J_to_eV * eV_to_qeV                # Combined conversion from Joules to qeV

    # --- Step 1: Unit Conversions and Reduced Mass (mu) ---
    r0 = r0_pm * pm_to_m
    m_H = m_H_amu * amu_to_kg
    m_Br = m_Br_amu * amu_to_kg
    mu = (m_H * m_Br) / (m_H + m_Br)

    # --- Step 2: Calculate Intermediate Constants (in SI units) ---
    # Moment of Inertia (I)
    I = mu * r0**2

    # Rotational Constant (Be) in Joules
    Be_J = h**2 / (8 * (pi**2) * I)

    # Vibrational Frequency (nu0) in Hz
    nu0 = (1 / (2 * pi)) * math.sqrt(k_force / mu)

    # Centrifugal Distortion Constant (De) in Joules
    De_J = (4 * Be_J**3) / ((h * nu0)**2)

    print("--- Calculated Molecular Constants ---")
    print(f"Reduced Mass (μ): {mu:.4e} kg")
    print(f"Moment of Inertia (I): {I:.4e} kg m^2")
    print(f"Centrifugal Distortion Constant (De): {De_J:.4e} J\n")
    print("--------------------------------------\n")
    
    # --- Step 3: Calculate Energy Shifts for Transitions ---
    print("The energy shift (ΔE) is given by: ΔE = -De * [J_final²(J_final+1)² - J_initial²(J_initial+1)²]\n")

    # Transition 1: J = 0 to J = 1
    J_i_1, J_f_1 = 0, 1
    term_1 = J_f_1**2 * (J_f_1 + 1)**2 - J_i_1**2 * (J_i_1 + 1)**2
    delta_E1_J = -De_J * term_1
    delta_E1_qeV = delta_E1_J * J_to_qeV
    
    print("1. For the rotational transition from J = 0 to J = 1:")
    print(f"   ΔE = -({De_J:.4e} J) * [{J_f_1}²*({J_f_1}+1)² - {J_i_1}²*({J_i_1}+1)²]")
    print(f"   ΔE = -({De_J:.4e} J) * [{int(term_1)}]")
    print(f"   ΔE = {delta_E1_J:.4e} J")
    print(f"   In quecto-electronvolts, this shift is: {delta_E1_qeV:.2e} qeV\n")

    # Transition 2: J = 1 to J = 2
    J_i_2, J_f_2 = 1, 2
    term_2 = J_f_2**2 * (J_f_2 + 1)**2 - J_i_2**2 * (J_i_2 + 1)**2
    delta_E2_J = -De_J * term_2
    delta_E2_qeV = delta_E2_J * J_to_qeV
    
    print("2. For the rotational transition from J = 1 to J = 2:")
    print(f"   ΔE = -({De_J:.4e} J) * [{J_f_2}²*({J_f_2}+1)² - {J_i_2}²*({J_i_2}+1)²]")
    print(f"   ΔE = -({De_J:.4e} J) * [{int(term_2)}]")
    print(f"   ΔE = {delta_E2_J:.4e} J")
    print(f"   In quecto-electronvolts, this shift is: {delta_E2_qeV:.2e} qeV\n")

# Run the calculation and print the results
calculate_energy_shifts()