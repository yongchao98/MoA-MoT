import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- 1. Constants and Input Parameters ---
    # Given data
    H_MASS_AMU = 1.008
    BR_MASS_AMU = 79.904
    R0_PM = 141.4
    K_N_M = 400.0

    # Physical constants
    AMU_TO_KG = 1.660539e-27  # kg/amu
    PM_TO_M = 1e-12           # m/pm
    H_PLANCK = 6.62607015e-34 # J*s
    H_BAR = H_PLANCK / (2 * math.pi) # J*s (Reduced Planck constant)
    E_CHARGE = 1.602176634e-19 # C (Elementary charge for eV conversion)

    # --- 2. Step-by-step Calculations ---
    # a. Reduced mass (mu) in kg
    m_h_kg = H_MASS_AMU * AMU_TO_KG
    m_br_kg = BR_MASS_AMU * AMU_TO_KG
    mu = (m_h_kg * m_br_kg) / (m_h_kg + m_br_kg)

    # b. Bond length (r0) in m
    r0_m = R0_PM * PM_TO_M

    # c. Moment of inertia (I) in kg*m^2
    I = mu * r0_m**2

    # d. Rotational constant (B) in Joules
    # B = ħ² / (2I)
    B_joules = (H_BAR**2) / (2 * I)

    # e. Angular vibrational frequency (omega) in rad/s
    # ω = sqrt(k / μ)
    omega_rad_s = math.sqrt(K_N_M / mu)

    # f. Centrifugal distortion constant (D) in Joules
    # D = 4B³ / (ħ²ω²)
    D_joules = (4 * B_joules**3) / (H_BAR**2 * omega_rad_s**2)

    # --- 3. Calculate Energy Shifts for Transitions ---
    # The energy shift for a transition J -> J+1 is ΔE = -4 * D * (J+1)³

    # For transition J=0 -> J=1, the initial state is J=0
    j_initial_1 = 0
    delta_E_joules_1 = -4 * D_joules * (j_initial_1 + 1)**3

    # For transition J=1 -> J=2, the initial state is J=1
    j_initial_2 = 1
    delta_E_joules_2 = -4 * D_joules * (j_initial_2 + 1)**3

    # --- 4. Convert to quecto-electronvolts (qeV) ---
    # 1 qeV = 10⁻³⁰ eV
    # 1 eV = E_CHARGE Joules
    # Conversion factor from Joules to qeV
    joules_to_qev = 1.0 / (1e-30 * E_CHARGE)

    delta_E_qev_1 = delta_E_joules_1 * joules_to_qev
    delta_E_qev_2 = delta_E_joules_2 * joules_to_qev

    # --- 5. Print the results ---
    print("Results of the calculation:\n")

    print("1. Energy shift for the rotational transition from J = 0 to J = 1:")
    print("   Equation: ΔE = -4 * D * (J+1)³")
    print(f"   Using D = {D_joules:.4e} J and J = {j_initial_1}:")
    print(f"   ΔE = -4 * {D_joules:.4e} J * ({j_initial_1}+1)³ = {delta_E_joules_1:.4e} J")
    print(f"   Energy Shift = {delta_E_qev_1:.3f} qeV")
    
    print("\n" + "-"*50 + "\n")
    
    print("2. Energy shift for the rotational transition from J = 1 to J = 2:")
    print("   Equation: ΔE = -4 * D * (J+1)³")
    print(f"   Using D = {D_joules:.4e} J and J = {j_initial_2}:")
    print(f"   ΔE = -4 * {D_joules:.4e} J * ({j_initial_2}+1)³ = {delta_E_joules_2:.4e} J")
    print(f"   Energy Shift = {delta_E_qev_2:.3f} qeV")


if __name__ == '__main__':
    calculate_energy_shifts()