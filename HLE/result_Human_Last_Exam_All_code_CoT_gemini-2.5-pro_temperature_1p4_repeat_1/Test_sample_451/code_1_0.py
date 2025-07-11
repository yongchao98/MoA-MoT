import math

def calculate_centrifugal_distortion_shift():
    """
    Calculates the energy shifts due to the centrifugal distortion effect
    for H-Br rotational transitions.
    """
    # --- 1. Define Constants and Given Values ---
    # Given parameters
    r0_pm = 141.4  # Bond length in picometers
    k_Nm = 400.0   # Force constant in N/m
    m_H_amu = 1.008  # Atomic mass of Hydrogen in amu
    m_Br_amu = 79.904 # Atomic mass of Bromine in amu

    # Physical constants
    h = 6.62607015e-34  # Planck's constant in J*s
    amu_to_kg = 1.66053906660e-27 # Conversion factor from amu to kg
    pm_to_m = 1e-12  # Conversion factor from pm to m
    eV_to_J = 1.602176634e-19  # Conversion factor from eV to J
    qeV_to_eV = 1e-30 # Conversion factor from qeV to eV

    # --- 2. Perform Intermediate Calculations in SI units ---
    # Convert given values to SI units
    r0_m = r0_pm * pm_to_m
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg

    # Calculate reduced mass (mu)
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Calculate moment of inertia (I)
    I_kg_m2 = mu_kg * r0_m**2

    # Calculate vibrational frequency (nu_e) in Hz
    # nu_e = (1 / 2*pi) * sqrt(k / mu)
    nu_e_Hz = (1 / (2 * math.pi)) * math.sqrt(k_Nm / mu_kg)

    # Calculate rotational constant (B) in Hz
    # B = h / (8 * pi^2 * I)
    B_Hz = h / (8 * math.pi**2 * I_kg_m2)

    # Calculate centrifugal distortion constant (D) in Hz
    # D = 4 * B^3 / nu_e^2
    D_Hz = (4 * B_Hz**3) / (nu_e_Hz**2)
    
    print("Analysis for H-Br molecule:")
    print(f"  Reduced mass (μ): {mu_kg:.5e} kg")
    print(f"  Moment of Inertia (I): {I_kg_m2:.5e} kg*m^2")
    print(f"  Vibrational Frequency (νe): {nu_e_Hz:.5e} Hz")
    print(f"  Rotational Constant (B): {B_Hz:.5e} Hz")
    print(f"  Centrifugal Distortion Constant (D): {D_Hz:.5e} Hz")
    print("-" * 50)

    # --- 3. Calculate Energy Shift for J=0 -> J=1 ---
    J_i_1, J_f_1 = 0, 1
    term_1 = J_f_1**2 * (J_f_1 + 1)**2 - J_i_1**2 * (J_i_1 + 1)**2
    
    delta_E_J_1 = -h * D_Hz * term_1
    delta_E_qeV_1 = delta_E_J_1 / eV_to_J / qeV_to_eV
    
    print("Calculation for transition J=0 -> J=1:")
    print("ΔE = -h * D * [J_f²(J_f+1)² - J_i²(J_i+1)²]")
    print(f"ΔE = -({h:.6e} J·s) * ({D_Hz:.6e} s⁻¹) * [{J_f_1}²({J_f_1}+1)² - {J_i_1}²({J_i_1}+1)²]")
    print(f"ΔE = -({h:.6e}) * ({D_Hz:.6e}) * [{int(term_1)}]")
    print(f"ΔE = {delta_E_J_1:.5e} J")
    print(f"ΔE = {delta_E_qeV_1:.5e} qeV")
    print("-" * 50)
    
    # --- 4. Calculate Energy Shift for J=1 -> J=2 ---
    J_i_2, J_f_2 = 1, 2
    term_2 = J_f_2**2 * (J_f_2 + 1)**2 - J_i_2**2 * (J_i_2 + 1)**2
    
    delta_E_J_2 = -h * D_Hz * term_2
    delta_E_qeV_2 = delta_E_J_2 / eV_to_J / qeV_to_eV

    print("Calculation for transition J=1 -> J=2:")
    print("ΔE = -h * D * [J_f²(J_f+1)² - J_i²(J_i+1)²]")
    print(f"ΔE = -({h:.6e} J·s) * ({D_Hz:.6e} s⁻¹) * [{J_f_2}²({J_f_2}+1)² - {J_i_2}²({J_i_2}+1)²]")
    print(f"ΔE = -({h:.6e}) * ({D_Hz:.6e}) * [{int(term_2)}]")
    print(f"ΔE = {delta_E_J_2:.5e} J")
    print(f"ΔE = {delta_E_qeV_2:.5e} qeV")

    return delta_E_qeV_1, delta_E_qeV_2

if __name__ == '__main__':
    ans1, ans2 = calculate_centrifugal_distortion_shift()
    # The final answer is provided by the detailed printout above.
    # The return values can be captured if this function is used elsewhere.

<<<[-1.76743e+23, -1.41394e+24]>>>