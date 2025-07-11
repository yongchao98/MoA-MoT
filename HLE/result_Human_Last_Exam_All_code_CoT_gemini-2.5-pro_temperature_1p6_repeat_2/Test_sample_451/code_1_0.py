import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br
    for J=0->1 and J=1->2 transitions.
    """
    # --- Step 1: Define Constants and Parameters ---
    # Given parameters
    m_H_amu = 1.008      # amu
    m_Br_amu = 79.904    # amu
    r0_pm = 141.4        # pm
    k_Nm = 400.0         # N/m

    # Physical constants
    amu_to_kg = 1.66053906660e-27  # kg
    pm_to_m = 1e-12                # m
    h = 6.62607015e-34             # J*s (Planck constant)
    eV_to_J = 1.602176634e-19      # J
    qeV_to_J = 1e-30 * eV_to_J     # J, for 1 quecto-electronvolt

    # Convert inputs to SI units
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    r0_m = r0_pm * pm_to_m

    # --- Step 2: Calculate Intermediate Molecular Properties ---
    # Reduced Mass (mu)
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Moment of Inertia (I)
    I_kg_m2 = mu_kg * r0_m**2

    # Rotational Constant (B_Hz) in Hz
    B_Hz = h / (8 * math.pi**2 * I_kg_m2)

    # Vibrational Frequency (nu_Hz) in Hz
    nu_Hz = (1 / (2 * math.pi)) * math.sqrt(k_Nm / mu_kg)

    # Centrifugal Distortion Constant (D_Hz) in Hz
    D_Hz = 4 * (B_Hz**3) / (nu_Hz**2)

    # --- Step 3 & 4: Calculate Energy Shifts and Convert Units ---

    # --- Transition 1: J = 0 -> J = 1 ---
    J_initial_1 = 0
    J_final_1 = 1
    # The term in the shift equation is J_final^2*(J_final+1)^2 - J_initial^2*(J_initial+1)^2
    term_1 = J_final_1**2 * (J_final_1 + 1)**2 - J_initial_1**2 * (J_initial_1 + 1)**2
    delta_E_J_1 = -h * D_Hz * term_1
    delta_E_qeV_1 = delta_E_J_1 / qeV_to_J

    # --- Transition 2: J = 1 -> J = 2 ---
    J_initial_2 = 1
    J_final_2 = 2
    term_2 = J_final_2**2 * (J_final_2 + 1)**2 - J_initial_2**2 * (J_initial_2 + 1)**2
    delta_E_J_2 = -h * D_Hz * term_2
    delta_E_qeV_2 = delta_E_J_2 / qeV_to_J

    # --- Step 5: Output Results ---
    print("--- Calculation for the J=0 to J=1 transition ---")
    print(f"The energy shift is given by the formula: \u0394E = -h \u00B7 D \u00B7 [ J_f\u00b2(J_f+1)\u00b2 - J_i\u00b2(J_i+1)\u00b2 ]")
    print(f"For J=0 \u2192 J=1, this becomes \u0394E = -h \u00B7 D \u00B7 [ {J_final_1**2 * (J_final_1 + 1)**2} - {J_initial_1**2 * (J_initial_1 + 1)**2} ] = -h \u00B7 D \u00B7 {term_1}")
    print(f"With h = {h:.5e} J\u00B7s and the calculated D = {D_Hz:.5e} Hz:")
    print(f"The energy shift \u0394E = -({h:.5e}) \u00B7 ({D_Hz:.5e}) \u00B7 {term_1:.1f} = {delta_E_J_1:.5e} J")
    print(f"In quecto-electronvolts, the shift is {delta_E_qeV_1:.4f} qeV.\n")

    print("--- Calculation for the J=1 to J=2 transition ---")
    print(f"The energy shift is given by the formula: \u0394E = -h \u00B7 D \u00B7 [ J_f\u00b2(J_f+1)\u00b2 - J_i\u00b2(J_i+1)\u00b2 ]")
    print(f"For J=1 \u2192 J=2, this becomes \u0394E = -h \u00B7 D \u00B7 [ {J_final_2**2 * (J_final_2 + 1)**2} - {J_initial_2**2 * (J_initial_2 + 1)**2} ] = -h \u00B7 D \u00B7 {term_2}")
    print(f"With h = {h:.5e} J\u00B7s and the calculated D = {D_Hz:.5e} Hz:")
    print(f"The energy shift \u0394E = -({h:.5e}) \u00B7 ({D_Hz:.5e}) \u00B7 {term_2:.1f} = {delta_E_J_2:.5e} J")
    print(f"In quecto-electronvolts, the shift is {delta_E_qeV_2:.4f} qeV.")
    
    # Returning final values for the <<<answer>>> tag
    return delta_E_qeV_1, delta_E_qeV_2

# Execute the function and capture the results
shift_1, shift_2 = calculate_centrifugal_distortion_shifts()

# The final answer is wrapped in <<<>>>
# print(f"\n<<<{shift_1:.4f}, {shift_2:.4f}>>>")