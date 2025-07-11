import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define Physical Constants
    H_PLANCK = 6.62607015e-34  # Planck's constant in J·s
    AMU_TO_KG = 1.66053906660e-27  # Atomic mass unit in kg
    E_CHARGE = 1.602176634e-19   # Elementary charge in C
    
    # Given molecular parameters
    r0 = 141.4e-12  # Bond length in meters
    k_force = 400.0  # Force constant in N/m
    m_H_amu = 1.008  # Mass of Hydrogen in amu
    m_Br_amu = 79.904 # Mass of Bromine in amu

    # 2. Calculate Molecular Properties
    # Convert masses to kg
    m_H_kg = m_H_amu * AMU_TO_KG
    m_Br_kg = m_Br_amu * AMU_TO_KG
    
    # Calculate reduced mass (mu)
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)
    
    # Calculate moment of inertia (I)
    I = mu * r0**2
    
    # 3. Determine Spectroscopic Constants
    # Calculate rotational constant (B) in Joules
    B_joules = H_PLANCK**2 / (8 * math.pi**2 * I)
    
    # Calculate vibrational frequency (nu) in Hz
    nu = (1 / (2 * math.pi)) * math.sqrt(k_force / mu)
    
    # Calculate vibrational energy quantum (h*nu) in Joules
    h_nu = H_PLANCK * nu
    
    # Calculate centrifugal distortion constant (D) in Joules
    # Using the approximation D ≈ 4B³ / (hν)²
    D_joules = 4 * B_joules**3 / h_nu**2
    
    # 4. Calculate Energy Shifts in Joules
    # For J=0 -> J=1 transition, the shift is -4D
    delta_E_0_to_1_J = -4 * D_joules
    
    # For J=1 -> J=2 transition, the shift is -32D
    delta_E_1_to_2_J = -32 * D_joules
    
    # 5. Convert Units and Output
    # Conversion factor from Joules to quecto-electronvolts (qeV)
    # 1 qeV = 1e-30 eV = 1e-30 * e_charge J
    J_TO_QEV = 1.0 / (E_CHARGE * 1e-30)
    
    delta_E_0_to_1_qeV = delta_E_0_to_1_J * J_TO_QEV
    delta_E_1_to_2_qeV = delta_E_1_to_2_J * J_TO_QEV
    
    print("Calculation Steps and Intermediate Values:")
    print(f"Reduced mass (μ) = {mu:.4e} kg")
    print(f"Moment of inertia (I) = {I:.4e} kg·m²")
    print(f"Rotational constant (B) = {B_joules:.4e} J")
    print(f"Vibrational frequency (ν) = {nu:.4e} Hz")
    print(f"Centrifugal distortion constant (D) = {D_joules:.4e} J")
    print("-" * 30)
    print("Final Energy Shifts:")
    print(f"The energy shift for the J=0 to J=1 transition is {delta_E_0_to_1_qeV:.4e} qeV.")
    print(f"The energy shift for the J=1 to J=2 transition is {delta_E_1_to_2_qeV:.4e} qeV.")

calculate_centrifugal_distortion_shifts()