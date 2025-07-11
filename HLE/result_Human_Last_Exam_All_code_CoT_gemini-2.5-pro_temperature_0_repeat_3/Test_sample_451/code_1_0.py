import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br
    for rotational transitions J=0->1 and J=1->2.
    """
    # 1. Constants and given values
    m_H_amu = 1.008         # amu
    m_Br_amu = 79.904       # amu
    r0_pm = 141.4           # pm
    k_Nm = 400              # N/m

    # Physical constants in SI units
    AMU_to_KG = 1.660539e-27  # kg/amu
    H_BAR = 1.0545718e-34     # J*s
    E_CHARGE = 1.6021766e-19  # J/eV
    PM_to_M = 1e-12           # m/pm
    QEV_per_EV = 1e30         # qeV/eV

    # Convert inputs to SI units
    m_H_kg = m_H_amu * AMU_to_KG
    m_Br_kg = m_Br_amu * AMU_to_KG
    r0_m = r0_pm * PM_to_M

    # 2. Calculate reduced mass (mu)
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # 3. Calculate moment of inertia (I)
    I = mu_kg * r0_m**2

    # 4. Calculate rotational constant (B) in Joules
    B_J = (H_BAR**2) / (2 * I)

    # 5. Calculate angular vibrational frequency (omega_e)
    omega_e = math.sqrt(k_Nm / mu_kg)

    # 6. Calculate centrifugal distortion constant (D) in Joules
    # Formula: D = 4 * B^3 / (hbar * omega_e)^2
    D_J = (4 * B_J**3) / ((H_BAR * omega_e)**2)

    # --- Calculations for Transitions ---

    # 7. Transition 1: J = 0 to J = 1
    # The energy shift is ΔE = -D * [J_f^2*(J_f+1)^2 - J_i^2*(J_i+1)^2]
    # For J=0->1, ΔE = -D * [1^2*(2)^2 - 0] = -4*D
    delta_E1_J = -4 * D_J

    # 8. Convert to qeV
    delta_E1_eV = delta_E1_J / E_CHARGE
    delta_E1_qeV = delta_E1_eV * QEV_per_EV

    print("Energy shift for J = 0 to J = 1 transition:")
    print(f"The centrifugal distortion constant is D = {D_J:.4e} J.")
    print(f"The energy shift equation is ΔE = -4 * D.")
    print(f"ΔE = -4 * {D_J:.4e} J = {delta_E1_J:.4e} J.")
    print(f"In electronvolts, ΔE = {delta_E1_J:.4e} J / {E_CHARGE:.4e} J/eV = {delta_E1_eV:.4e} eV.")
    print(f"The final energy shift is {delta_E1_eV:.4e} eV * {QEV_per_EV:.0e} qeV/eV = {delta_E1_qeV:.4e} qeV.")
    print("-" * 50)

    # 7. Transition 2: J = 1 to J = 2
    # For J=1->2, ΔE = -D * [2^2*(3)^2 - 1^2*(2)^2] = -D * [36 - 4] = -32*D
    delta_E2_J = -32 * D_J

    # 8. Convert to qeV
    delta_E2_eV = delta_E2_J / E_CHARGE
    delta_E2_qeV = delta_E2_eV * QEV_per_EV

    print("Energy shift for J = 1 to J = 2 transition:")
    print(f"The centrifugal distortion constant is D = {D_J:.4e} J.")
    print(f"The energy shift equation is ΔE = -32 * D.")
    print(f"ΔE = -32 * {D_J:.4e} J = {delta_E2_J:.4e} J.")
    print(f"In electronvolts, ΔE = {delta_E2_J:.4e} J / {E_CHARGE:.4e} J/eV = {delta_E2_eV:.4e} eV.")
    print(f"The final energy shift is {delta_E2_eV:.4e} eV * {QEV_per_EV:.0e} qeV/eV = {delta_E2_qeV:.4e} qeV.")
    
    # Format final answer for submission
    global final_answer
    final_answer = f"For J=0->1: {delta_E1_qeV:.3e} qeV, For J=1->2: {delta_E2_qeV:.3e} qeV"


if __name__ == '__main__':
    calculate_centrifugal_distortion_shifts()
    # The final answer is captured in a global variable to be appended at the end.
    # print(f"\n<<<{final_answer}>>>")