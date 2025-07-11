import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # Step 1: Define Constants
    # Given data
    m_H_amu = 1.008       # atomic mass units (amu)
    m_Br_amu = 79.904     # amu
    r0_pm = 141.4         # picometers (pm)
    k_Nm = 400            # Newtons per meter (N/m)

    # Physical constants in SI units
    amu_to_kg = 1.66053906660e-27  # kg/amu
    pm_to_m = 1e-12                # m/pm
    hbar = 1.054571817e-34         # J·s
    e_charge = 1.602176634e-19     # Coulombs (for J to eV conversion)

    # Step 2: Calculate Reduced Mass (mu) in kg
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Convert bond length to meters
    r0_m = r0_pm * pm_to_m

    # Step 3: Calculate Centrifugal Distortion Constant (D) in Joules
    # D = hbar^4 / (2 * mu^2 * r0^6 * k)
    D_J = (hbar**4) / (2 * mu_kg**2 * r0_m**6 * k_Nm)

    # Step 4: Calculate Energy Shift for each transition in Joules
    # For J=0 to J=1 transition, the shift is -4D
    delta_E1_J = -4 * D_J

    # For J=1 to J=2 transition, the shift is -32D
    delta_E2_J = -32 * D_J

    # Step 5: Convert from Joules to quecto-electronvolts (qeV)
    # 1 qeV = 1e-30 eV
    # To convert J to qeV: (Value in J / e_charge) * 1e30
    joules_to_qev_factor = 1e30 / e_charge
    
    delta_E1_qeV = delta_E1_J * joules_to_qev_factor
    delta_E2_qeV = delta_E2_J * joules_to_qev_factor

    # Print the results
    print("This script calculates the energy shifts due to centrifugal distortion for H-Br.")
    print("-" * 70)
    print(f"Calculated Centrifugal Distortion Constant (D) = {D_J:.4e} J")
    print("-" * 70)
    
    # Final Equation for Transition 1: Energy Shift = -4 * D
    print("1. For the J=0 to J=1 transition:")
    print(f"   Energy Shift = -4 * {D_J:.4e} J = {delta_E1_J:.4e} J")
    print(f"   In quecto-electronvolts, this is: {delta_E1_qeV:.3e} qeV\n")

    # Final Equation for Transition 2: Energy Shift = -32 * D
    print("2. For the J=1 to J=2 transition:")
    print(f"   Energy Shift = -32 * {D_J:.4e} J = {delta_E2_J:.4e} J")
    print(f"   In quecto-electronvolts, this is: {delta_E2_qeV:.3e} qeV")

calculate_energy_shifts()