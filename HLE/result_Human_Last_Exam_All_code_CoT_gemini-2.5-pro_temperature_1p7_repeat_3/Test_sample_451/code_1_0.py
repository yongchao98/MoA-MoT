import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br.
    """
    # 1. Define constants and parameters
    # Physical constants
    HBAR = 1.054571817e-34  # Reduced Planck constant (J*s)
    AMU_TO_KG = 1.66053906660e-27  # Atomic mass unit to kg conversion
    PM_TO_M = 1e-12  # Picometer to meter conversion
    E_CHARGE = 1.602176634e-19 # Elementary charge in Coulombs
    
    # Molecule-specific parameters for H-Br
    m_H_amu = 1.008         # Atomic mass of Hydrogen (amu)
    m_Br_amu = 79.904       # Atomic mass of Bromine (amu)
    r0_pm = 141.4           # Bond length (pm)
    k_force = 400.0         # Force constant (N/m)

    # Unit conversions
    r0 = r0_pm * PM_TO_M
    m_H = m_H_amu * AMU_TO_KG
    m_Br = m_Br_amu * AMU_TO_KG
    
    # 2. Calculate molecular properties
    # Reduced mass (mu)
    mu = (m_H * m_Br) / (m_H + m_Br)
    
    # Moment of inertia (I)
    I = mu * r0**2
    
    # 3. Calculate spectroscopic constants in SI units (Joules)
    # Rotational constant (B) in Joules
    B_J = (HBAR**2) / (2 * I)
    
    # Fundamental vibrational energy (hv0) in Joules
    # Angular frequency omega = sqrt(k/mu)
    # hv0 = hbar * omega
    h_nu_J = HBAR * math.sqrt(k_force / mu)
    
    # Centrifugal distortion constant (D) in Joules
    D_J = (4 * B_J**3) / (h_nu_J**2)

    print("--- Intermediate Calculations (in SI Units) ---")
    print(f"Reduced mass (μ): {mu:.4e} kg")
    print(f"Moment of inertia (I): {I:.4e} kg m^2")
    print(f"Rotational constant (B): {B_J:.4e} J")
    print(f"Vibrational energy quantum (hν₀): {h_nu_J:.4e} J")
    print(f"Centrifugal distortion constant (D): {D_J:.4e} J")
    print("-" * 40)
    
    # 4. & 5. Calculate energy shifts and convert units
    # Conversion from Joules to qeV
    # 1 qeV = 1e-30 eV = 1e-30 * e_charge J
    J_PER_QEV = 1e-30 * E_CHARGE
    
    print("\n--- Energy Shift Calculations ---")

    # Transition 1: J = 0 -> J = 1
    J1_initial = 0
    shift_J_0_1 = -4 * D_J * (J1_initial + 1)**3
    shift_qeV_0_1 = shift_J_0_1 / J_PER_QEV
    
    print("\n1. For the transition J = 0 to J = 1:")
    print(f"   ΔE_shift = -4 * D * (J+1)^3")
    print(f"   ΔE_shift = -4 * ({D_J:.4e} J) * ({J1_initial}+1)^3")
    print(f"   ΔE_shift = {shift_J_0_1:.4e} J")
    print(f"   Converting to qeV: {shift_J_0_1:.4e} J / ({J_PER_QEV:.4e} J/qeV)")
    print(f"   Energy shift = {shift_qeV_0_1:.4e} qeV")
    
    # Transition 2: J = 1 -> J = 2
    J2_initial = 1
    shift_J_1_2 = -4 * D_J * (J2_initial + 1)**3
    shift_qeV_1_2 = shift_J_1_2 / J_PER_QEV
    
    print("\n2. For the transition J = 1 to J = 2:")
    print(f"   ΔE_shift = -4 * D * (J+1)^3")
    print(f"   ΔE_shift = -4 * ({D_J:.4e} J) * ({J2_initial}+1)^3")
    print(f"   ΔE_shift = {shift_J_1_2:.4e} J")
    print(f"   Converting to qeV: {shift_J_1_2:.4e} J / ({J_PER_QEV:.4e} J/qeV)")
    print(f"   Energy shift = {shift_qeV_1_2:.4e} qeV")

calculate_centrifugal_distortion_shifts()
<<<The energy shift for the J=0->1 transition is -1.7672e+23 qeV and for the J=1->2 transition is -1.4137e+24 qeV.>>>