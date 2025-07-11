import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # Step 1: Define Constants and Inputs
    # Physical constants
    h_cgs = 6.62607015e-27      # Planck's constant in erg*s
    h_si = 6.62607015e-34       # Planck's constant in J*s
    c_cgs = 2.99792458e10       # Speed of light in cm/s
    c_si = 2.99792458e8         # Speed of light in m/s
    amu_to_kg = 1.66053906660e-27 # amu to kg conversion
    e_charge = 1.602176634e-19   # Elementary charge in Coulombs

    # Given data
    r0_pm = 141.4               # Bond length in pm
    k_si = 400.0                # Force constant in N/m
    m_H_amu = 1.008             # Mass of Hydrogen in amu
    m_Br_amu = 79.904           # Mass of Bromine in amu

    # Convert inputs to SI units
    r0_m = r0_pm * 1e-12

    # Step 2: Calculate Molecular Properties
    # Reduced mass in kg
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Moment of inertia in kg*m^2
    I_si = mu_kg * r0_m**2

    # Step 3: Calculate Spectroscopic Constants in cm^-1
    # Vibrational frequency in cm^-1
    # omega_e (rad/s) = sqrt(k/mu)
    # omega_e (cm^-1) = omega_e(rad/s) / (2 * pi * c_cm/s)
    omega_e_cm = (1 / (2 * math.pi * c_cgs)) * math.sqrt(k_si / mu_kg)

    # Rotational constant in cm^-1
    # B_e (J) = h_si^2 / (8 * pi^2 * I_si)
    # B_e (cm^-1) = B_e(J) / (h_cgs * c_cgs) = h_si / (8 * pi^2 * c_cgs * I_si)
    B_e_cm = h_si / (8 * math.pi**2 * c_cgs * I_si)

    # Centrifugal distortion constant in cm^-1
    D_e_cm = (4 * (B_e_cm**3)) / (omega_e_cm**2)

    # Step 4: Calculate Energy Shifts in Joules
    # Convert D_e from cm^-1 to Joules
    # Energy (J) = Energy (cm^-1) * h (J.s) * c (m/s) * 100 (cm/m)
    D_J = D_e_cm * h_si * c_si * 100

    # Calculate magnitude of energy shift for J -> J+1: |ΔE| = 4 * D * (J+1)^3
    # Transition 1: J = 0 -> J = 1
    J_initial_1 = 0
    delta_E_J_1 = 4 * D_J * (J_initial_1 + 1)**3

    # Transition 2: J = 1 -> J = 2
    J_initial_2 = 1
    delta_E_J_2 = 4 * D_J * (J_initial_2 + 1)**3

    # Step 5: Convert to Final Units (qeV) and Print
    # 1 qeV = 1e-30 eV = 1e-30 * e_charge Joules
    joules_to_qev = 1.0 / (1e-30 * e_charge)

    delta_E_qeV_1 = delta_E_J_1 * joules_to_qev
    delta_E_qeV_2 = delta_E_J_2 * joules_to_qev

    # Print the results
    print(f"The energy shift due to centrifugal distortion is calculated using the formula |ΔE_shift| = 4 * D * (J+1)³.")
    print(f"The calculated centrifugal distortion constant, D, is {D_J:.4e} Joules.")
    print("\n--- Results ---")
    print(f"For the transition from J = {J_initial_1} to J = {J_initial_1 + 1}:")
    print(f"The energy shift is 4 * ({D_J:.4e} J) * ({J_initial_1} + 1)³ = {delta_E_qeV_1:.3f} qeV")
    
    print(f"\nFor the transition from J = {J_initial_2} to J = {J_initial_2 + 1}:")
    print(f"The energy shift is 4 * ({D_J:.4e} J) * ({J_initial_2} + 1)³ = {delta_E_qeV_2:.3f} qeV")

calculate_energy_shifts()