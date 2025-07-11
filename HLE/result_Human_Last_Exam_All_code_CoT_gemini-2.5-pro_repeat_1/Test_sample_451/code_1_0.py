import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- 1. Define Constants and Input Parameters ---
    # Given parameters
    m_H_amu = 1.008
    m_Br_amu = 79.904
    r0_pm = 141.4
    k = 400.0  # N/m

    # Physical constants
    amu_to_kg = 1.660539e-27
    pm_to_m = 1e-12
    h = 6.62607015e-34       # Planck's constant in J*s
    eV_to_J = 1.602176634e-19 # J per eV
    qeV_to_eV = 1e-30         # eV per qeV

    # --- 2. Calculate Molecular Properties and Spectroscopic Constants ---
    # Convert masses to kg
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg

    # Calculate reduced mass (μ) in kg
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Calculate classical vibrational frequency (ω) in Hz
    omega = (1 / (2 * math.pi)) * math.sqrt(k / mu)

    # Convert bond length to meters and calculate moment of inertia (I) in kg*m^2
    r0_m = r0_pm * pm_to_m
    I = mu * r0_m**2

    # Calculate rotational constant (B) in Hz
    B = h / (8 * math.pi**2 * I)

    # Calculate centrifugal distortion constant (D) in Hz
    D = 4 * B**3 / omega**2
    
    # Conversion factor from Joules to quecto-electronvolts
    J_to_qeV = 1 / (eV_to_J * qeV_to_eV)

    # --- 3. Calculate Energy Shift for J=0 -> J=1 ---
    print("--- Transition 1: J = 0 -> J = 1 ---")
    
    # Calculate energy shift in Joules
    # Formula: ΔE = -4 * h * D
    delta_E1_J = -4 * h * D
    
    print(f"The energy shift is calculated as ΔE = -4 * h * D")
    print(f"  h = {h:.4e} J*s")
    print(f"  D = {D:.4e} Hz")
    print(f"ΔE = -4 * ({h:.4e}) * ({D:.4e})")
    print(f"ΔE = {delta_E1_J:.4e} J")
    
    # Convert to qeV
    delta_E1_qeV = delta_E1_J * J_to_qeV
    print(f"\nThe energy shift for the J=0 -> J=1 transition is: {delta_E1_qeV:.4e} qeV\n")

    # --- 4. Calculate Energy Shift for J=1 -> J=2 ---
    print("--- Transition 2: J = 1 -> J = 2 ---")
    
    # Calculate energy shift in Joules
    # Formula: ΔE = -32 * h * D
    delta_E2_J = -32 * h * D
    
    print(f"The energy shift is calculated as ΔE = -32 * h * D")
    print(f"  h = {h:.4e} J*s")
    print(f"  D = {D:.4e} Hz")
    print(f"ΔE = -32 * ({h:.4e}) * ({D:.4e})")
    print(f"ΔE = {delta_E2_J:.4e} J")

    # Convert to qeV
    delta_E2_qeV = delta_E2_J * J_to_qeV
    print(f"\nThe energy shift for the J=1 -> J=2 transition is: {delta_E2_qeV:.4e} qeV")


if __name__ == '__main__':
    calculate_energy_shifts()