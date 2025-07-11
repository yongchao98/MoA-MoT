import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- 1. Constants and Given Data ---
    m_H_amu = 1.008
    m_Br_amu = 79.904
    r_0_pm = 141.4
    k_Nm = 400.0
    
    # Physical constants in SI units
    AMU_TO_KG = 1.66053906660e-27
    PM_TO_M = 1e-12
    H_JS = 6.62607015e-34
    E_CHARGE = 1.602176634e-19
    J_TO_EV = 1 / E_CHARGE
    EV_TO_QEV = 1e30
    J_TO_QEV = J_TO_EV * EV_TO_QEV

    print("--- Calculating Input Parameters in SI units ---")
    
    # --- 2. Intermediate Calculations in SI units ---
    
    # Reduced mass (mu)
    m_H_kg = m_H_amu * AMU_TO_KG
    m_Br_kg = m_Br_amu * AMU_TO_KG
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)
    print(f"Reduced mass (μ): {mu_kg:.4e} kg")

    # Bond length (r_0)
    r_0_m = r_0_pm * PM_TO_M
    print(f"Bond length (r_0): {r_0_m:.4e} m")
    
    # Moment of inertia (I)
    I_kg_m2 = mu_kg * r_0_m**2
    print(f"Moment of inertia (I): {I_kg_m2:.4e} kg m^2")
    
    # Rotational constant (B) in Hz
    B_hz = H_JS / (8 * math.pi**2 * I_kg_m2)
    print(f"Rotational constant (B): {B_hz:.4e} Hz")
    
    # Vibrational frequency (nu) in Hz
    nu_hz = (1 / (2 * math.pi)) * math.sqrt(k_Nm / mu_kg)
    print(f"Vibrational frequency (ν): {nu_hz:.4e} Hz")
    
    # Centrifugal distortion constant (D) in Hz
    D_hz = 4 * B_hz**3 / nu_hz**2
    print(f"Centrifugal distortion constant (D): {D_hz:.4e} Hz")
    
    # Energy term h*D in Joules
    hD_joules = H_JS * D_hz
    print(f"Energy term (h * D): {hD_joules:.4e} J\n")

    # --- 3. Calculate Energy Shift for J=0 to J=1 ---
    print("--- 1. Energy Shift for J = 0 -> 1 Transition ---")
    
    # The energy shift is ΔE = -4 * h * D
    shift_factor_1 = -4
    energy_shift_1_joules = shift_factor_1 * hD_joules
    energy_shift_1_qeV = energy_shift_1_joules * J_TO_QEV

    print(f"The formula for the energy shift is: ΔE = {shift_factor_1} * h * D")
    print("Plugging in the values:")
    print(f"ΔE(0->1) = {shift_factor_1} * ({hD_joules:.4e} J)")
    print(f"ΔE(0->1) = {energy_shift_1_joules:.4e} J")
    print(f"ΔE(0->1) = {energy_shift_1_qeV:.4e} qeV\n")

    # --- 4. Calculate Energy Shift for J=1 to J=2 ---
    print("--- 2. Energy Shift for J = 1 -> 2 Transition ---")
    
    # The energy shift is ΔE = -32 * h * D
    shift_factor_2 = -32
    energy_shift_2_joules = shift_factor_2 * hD_joules
    energy_shift_2_qeV = energy_shift_2_joules * J_TO_QEV

    print(f"The formula for the energy shift is: ΔE = {shift_factor_2} * h * D")
    print("Plugging in the values:")
    print(f"ΔE(1->2) = {shift_factor_2} * ({hD_joules:.4e} J)")
    print(f"ΔE(1->2) = {energy_shift_2_joules:.4e} J")
    print(f"ΔE(1->2) = {energy_shift_2_qeV:.4e} qeV")

calculate_energy_shifts()