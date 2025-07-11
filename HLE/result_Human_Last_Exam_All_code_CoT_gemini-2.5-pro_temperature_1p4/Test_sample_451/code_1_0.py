import math
import scipy.constants as const

def calculate_centrifugal_distortion_shift():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br
    for J=0->1 and J=1->2 transitions.
    """
    # 1. Define Constants and Inputs
    m_H_amu = 1.008  # amu
    m_Br_amu = 79.904 # amu
    r0 = 141.4e-12    # meters
    k = 400.0         # N/m
    
    # Physical constants from scipy
    m_u = const.m_u      # Atomic mass unit in kg
    h = const.h          # Planck constant in J*s
    e = const.e          # Elementary charge in C
    pi = math.pi

    # 2. Calculate Key Molecular Properties
    # Convert masses to kg
    m_H = m_H_amu * m_u
    m_Br = m_Br_amu * m_u
    
    # Reduced mass (mu)
    mu = (m_H * m_Br) / (m_H + m_Br)
    
    # Moment of inertia (I)
    I = mu * r0**2

    # 3. Determine Spectroscopic Constants
    # Rotational constant (B) in Hz
    B_Hz = h / (8 * pi**2 * I)
    
    # Vibrational frequency (nu0) in Hz
    nu0 = (1 / (2 * pi)) * math.sqrt(k / mu)
    
    # Centrifugal distortion constant (D) in Hz
    D_Hz = (4 * B_Hz**3) / (nu0**2)
    
    # The energy term h*D in Joules
    hD_J = h * D_Hz

    # 4. Calculate Energy Shifts in Joules
    # For J=0 -> J=1, the shift is -4 * h * D
    delta_E_0_1_J = -4 * hD_J
    
    # For J=1 -> J=2, the shift is -4 * h * D * (1+1)^3 = -32 * h * D
    delta_E_1_2_J = -32 * hD_J
    
    # 5. Convert Units and Output
    # Conversion factor from Joules to quecto-electronvolts (qeV)
    # 1 qeV = 1e-30 eV = 1e-30 * e Joules
    J_per_qeV = e * 1e-30
    
    delta_E_0_1_qeV = delta_E_0_1_J / J_per_qeV
    delta_E_1_2_qeV = delta_E_1_2_J / J_per_qeV

    print("--- Calculation of Energy Shifts for H-Br ---")
    
    print("\n1. For the rotational transition from J = 0 to J = 1:")
    print(f"Energy Shift = -4 * h * D")
    print(f"               = -4 * ({hD_J:.5e} J)")
    print(f"               = {delta_E_0_1_J:.5e} J")
    print(f"In qeV:        = ({delta_E_0_1_J:.5e} J) / ({J_per_qeV:.5e} J/qeV)")
    print(f"               = {delta_E_0_1_qeV:.4e} qeV")

    print("\n2. For the rotational transition from J = 1 to J = 2:")
    print(f"Energy Shift = -32 * h * D")
    print(f"               = -32 * ({hD_J:.5e} J)")
    print(f"               = {delta_E_1_2_J:.5e} J")
    print(f"In qeV:        = ({delta_E_1_2_J:.5e} J) / ({J_per_qeV:.5e} J/qeV)")
    print(f"               = {delta_E_1_2_qeV:.4e} qeV")

if __name__ == '__main__':
    calculate_centrifugal_distortion_shift()