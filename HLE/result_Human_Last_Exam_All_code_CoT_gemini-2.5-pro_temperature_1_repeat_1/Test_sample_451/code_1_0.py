import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define Constants and Input Values
    # Given parameters
    m_h_amu = 1.008       # Mass of Hydrogen (amu)
    m_br_amu = 79.904     # Mass of Bromine (amu)
    r0_pm = 141.4         # Bond length (pm)
    k_Nm = 400            # Force constant (N/m)

    # Physical constants
    AMU_TO_KG = 1.660539e-27  # Conversion factor for amu to kg
    HBAR = 1.0545718e-34      # Reduced Planck constant (J*s)
    E_CHARGE = 1.6021766e-19  # Elementary charge (C)
    PM_TO_M = 1e-12           # Conversion factor for pm to m
    J_TO_QEV = 1e30 / E_CHARGE # Conversion factor from Joules to quecto-electronvolts

    # 2. Unit Conversions
    m_h_kg = m_h_amu * AMU_TO_KG
    m_br_kg = m_br_amu * AMU_TO_KG
    r0_m = r0_pm * PM_TO_M

    # 3. Calculate Molecular Properties
    # Reduced Mass (mu)
    mu_kg = (m_h_kg * m_br_kg) / (m_h_kg + m_br_kg)

    # Moment of Inertia (I)
    I = mu_kg * r0_m**2

    # Rotational Constant (B) in Joules
    B_J = (HBAR**2) / (2 * I)

    # 4. Calculate Centrifugal Distortion Constant (D) in Joules
    # The formula is D = 4 * B^3 / (hbar * omega)^2, where omega^2 = k / mu
    # This simplifies to D = (4 * B^3 * mu) / (hbar^2 * k)
    D_J = (4 * B_J**3 * mu_kg) / (HBAR**2 * k_Nm)

    # 5. Calculate Energy Shifts for each transition
    # --- Transition 1: J = 0 to J = 1 ---
    J1 = 0
    # Formula for shift: delta_E = -4 * D * (J+1)^3
    delta_E1_J = -4 * D_J * (J1 + 1)**3
    
    # --- Transition 2: J = 1 to J = 2 ---
    J2 = 1
    # Formula for shift: delta_E = -4 * D * (J+1)^3
    delta_E2_J = -4 * D_J * (J2 + 1)**3
    
    # 6. Convert to qeV and Output Results
    delta_E1_qeV = delta_E1_J * J_TO_QEV
    delta_E2_qeV = delta_E2_J * J_TO_QEV

    print("--- H-Br Centrifugal Distortion Calculation ---")
    print(f"Calculated centrifugal distortion constant D = {D_J:.4e} J\n")

    # Output for J=0 to J=1 transition
    print("1. Energy shift for transition from J = 0 to J = 1:")
    print(f"   Formula: ΔE = -4 * D * (J+1)^3")
    print(f"   Calculation: ΔE = -4 * ({D_J:.4e} J) * ({J1}+1)^3 = {delta_E1_J:.4e} J")
    print(f"   Energy shift = {delta_E1_qeV:.4e} qeV\n")

    # Output for J=1 to J=2 transition
    print("2. Energy shift for transition from J = 1 to J = 2:")
    print(f"   Formula: ΔE = -4 * D * (J+1)^3")
    print(f"   Calculation: ΔE = -4 * ({D_J:.4e} J) * ({J2}+1)^3 = {delta_E2_J:.4e} J")
    print(f"   Energy shift = {delta_E2_qeV:.4e} qeV")

calculate_centrifugal_distortion_shifts()