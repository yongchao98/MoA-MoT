import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br
    for rotational transitions J=0->1 and J=1->2.
    """

    # 1. Define physical constants and input parameters
    H_PLANCK = 6.62607015e-34  # Planck's constant in J*s
    AMU_TO_KG = 1.660539e-27  # Atomic mass unit to kg conversion
    EV_TO_J = 1.602176634e-19  # Electronvolt to Joules conversion
    QEV_TO_J = EV_TO_J * 1e-30 # quecto-electronvolt to Joules conversion

    # Given parameters for H-Br
    r0 = 141.4e-12  # Bond length in meters (141.4 pm)
    k = 400.0       # Force constant in N/m
    m_H = 1.008     # Mass of Hydrogen in amu
    m_Br = 79.904   # Mass of Bromine in amu

    # 2. Calculate reduced mass (mu) in kg
    m_H_kg = m_H * AMU_TO_KG
    m_Br_kg = m_Br * AMU_TO_KG
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # 3. Calculate moment of inertia (I)
    I = mu * r0**2

    # 4. Calculate rotational constant (B) in Hz
    B = H_PLANCK / (8 * math.pi**2 * I)

    # 5. Calculate vibrational frequency (omega) in Hz
    omega = (1 / (2 * math.pi)) * math.sqrt(k / mu)

    # 6. Calculate centrifugal distortion constant (D) in Hz
    D = (4 * B**3) / (omega**2)

    # 7. Calculate energy shifts for each transition
    # The formula for the shift in transition energy is: ΔE = -4 * h * D * (J_lower + 1)^3

    # For transition J=0 to J=1:
    J_lower_1 = 0
    factor_1 = -4 * (J_lower_1 + 1)**3
    delta_E_joules_1 = factor_1 * H_PLANCK * D
    delta_E_qev_1 = delta_E_joules_1 / QEV_TO_J

    # For transition J=1 to J=2:
    J_lower_2 = 1
    factor_2 = -4 * (J_lower_2 + 1)**3
    delta_E_joules_2 = factor_2 * H_PLANCK * D
    delta_E_qev_2 = delta_E_joules_2 / QEV_TO_J

    # Print the detailed results
    print("This script calculates the energy shift due to centrifugal distortion for H-Br.")
    print("-" * 60)

    print("1. For the rotational transition from J = 0 to J = 1:")
    print(f"   The energy shift is given by the equation: ΔE = {factor_1} * h * D")
    print(f"   ΔE = {factor_1} * ({H_PLANCK:.4e} J·s) * ({D:.4e} Hz)")
    print(f"   The resulting energy shift is {delta_E_qev_1:.4e} qeV.")
    print("-" * 60)

    print("2. For the rotational transition from J = 1 to J = 2:")
    print(f"   The energy shift is given by the equation: ΔE = {factor_2} * h * D")
    print(f"   ΔE = {factor_2} * ({H_PLANCK:.4e} J·s) * ({D:.4e} Hz)")
    print(f"   The resulting energy shift is {delta_E_qev_2:.4e} qeV.")
    print("-" * 60)

if __name__ == '__main__':
    calculate_centrifugal_distortion_shifts()