import scipy.constants as const

def calculate_centrifugal_distortion_shift():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define Constants and Inputs
    m_H_amu = 1.008  # amu
    m_Br_amu = 79.904 # amu
    r0_pm = 141.4    # picometers
    k_force = 400.0  # N/m

    # Physical constants from scipy
    amu_to_kg = const.physical_constants['atomic mass unit-kilogram relationship'][0]
    hbar = const.hbar  # J·s
    eV_to_J = const.e  # J/eV

    # 2. Unit Conversions
    m_H = m_H_amu * amu_to_kg
    m_Br = m_Br_amu * amu_to_kg
    r0 = r0_pm * 1e-12  # meters

    # 3. Calculate Molecular Properties
    # Reduced mass (mu)
    mu = (m_H * m_Br) / (m_H + m_Br)

    # Moment of inertia (I)
    I = mu * r0**2

    # 4. Calculate Spectroscopic Constants (in Joules)
    # Rotational constant (B) in Joules
    B_joules = (hbar**2) / (2 * I)

    # Angular vibrational frequency (omega) in rad/s
    omega_rad_s = (k_force / mu)**0.5

    # Centrifugal distortion constant (D) in Joules
    # The energy quantum for the vibrator is hbar*omega
    D_joules = (4 * B_joules**3) / (hbar * omega_rad_s)**2

    # 5. Calculate Energy Shifts in Joules
    # For J=0 -> J=1 transition, shift is -4*D
    delta_E_0_1_J = -4 * D_joules
    
    # For J=1 -> J=2 transition, shift is -32*D
    delta_E_1_2_J = -32 * D_joules

    # 6. Convert to qeV and Print
    # Conversion factor from Joules to qeV (1 qeV = 1e-30 eV)
    qeV_to_J = 1e-30 * eV_to_J
    J_to_qeV = 1 / qeV_to_J
    
    # Convert D to qeV to show the components of the final equation
    D_qeV = D_joules * J_to_qeV

    delta_E_0_1_qeV = delta_E_0_1_J * J_to_qeV
    delta_E_1_2_qeV = delta_E_1_2_J * J_to_qeV

    print("Energy shift for the J=0 to J=1 transition:")
    print(f"ΔE(0->1) = -4 * D = -4 * ({D_qeV:.4e}) qeV = {delta_E_0_1_qeV:.4e} qeV")
    print("\nEnergy shift for the J=1 to J=2 transition:")
    print(f"ΔE(1->2) = -32 * D = -32 * ({D_qeV:.4e}) qeV = {delta_E_1_2_qeV:.4e} qeV")

calculate_centrifugal_distortion_shift()