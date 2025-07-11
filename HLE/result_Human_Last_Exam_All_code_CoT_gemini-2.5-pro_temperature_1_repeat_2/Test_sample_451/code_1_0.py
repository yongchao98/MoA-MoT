import scipy.constants as const

def calculate_hbr_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for rotational
    transitions in an H-Br molecule.
    """
    # --- Constants and Input Parameters ---
    # Given values
    r0_pm = 141.4      # Bond length in pm
    k = 400.0          # Force constant in N/m
    m_H_amu = 1.008    # Mass of Hydrogen in amu
    m_Br_amu = 79.904  # Mass of Bromine in amu

    # Physical constants
    h_bar = const.hbar      # Reduced Planck constant in J·s
    amu_to_kg = const.u     # Atomic mass unit in kg
    J_to_eV = 1.0 / const.e # Conversion factor from Joules to eV

    # Unit conversions
    r0 = r0_pm * 1e-12      # Convert pm to m
    m_H = m_H_amu * amu_to_kg
    m_Br = m_Br_amu * amu_to_kg
    eV_to_qeV = 1e30        # Conversion factor from eV to qeV

    # --- Calculations ---
    # 1. Calculate the reduced mass (mu)
    mu = (m_H * m_Br) / (m_H + m_Br)

    # 2. Calculate the centrifugal distortion constant (D) in Joules
    # Formula: D = h_bar^4 / (2 * mu^2 * r0^6 * k)
    D_joules = (h_bar**4) / (2 * mu**2 * r0**6 * k)

    # 3. Calculate energy shifts for the specified transitions
    # Formula for shift: Delta_E = -4 * D * (J_initial + 1)^3

    # Transition 1: J = 0 to J = 1
    J1_initial = 0
    shift1_joules = -4 * D_joules * (J1_initial + 1)**3
    shift1_qeV = shift1_joules * J_to_eV * eV_to_qeV

    # Transition 2: J = 1 to J = 2
    J2_initial = 1
    shift2_joules = -4 * D_joules * (J2_initial + 1)**3
    shift2_qeV = shift2_joules * J_to_eV * eV_to_qeV

    # --- Output Results ---
    print("Calculation of energy shifts for H-Br due to centrifugal distortion:")
    print("-" * 60)
    print("Intermediate calculated values:")
    print(f"  Reduced mass (μ) = {mu:.4e} kg")
    print(f"  Centrifugal distortion constant (D) = {D_joules:.4e} J")
    print("-" * 60)

    # Output for the J=0 to J=1 transition
    print("1. For the transition from J = 0 to J = 1:")
    print(f"   The energy shift is given by the equation: ΔE = -4 * D * (J+1)³")
    print(f"   Plugging in the values: ΔE = -4 * ({D_joules:.4e} J) * ({J1_initial}+1)³")
    print(f"   Resulting energy shift: ΔE = {shift1_joules:.4e} J")
    print(f"   In quecto-electronvolts, the energy shift is: {shift1_qeV:.4e} qeV")
    print()

    # Output for the J=1 to J=2 transition
    print("2. For the transition from J = 1 to J = 2:")
    print(f"   The energy shift is given by the equation: ΔE = -4 * D * (J+1)³")
    print(f"   Plugging in the values: ΔE = -4 * ({D_joules:.4e} J) * ({J2_initial}+1)³")
    print(f"   Resulting energy shift: ΔE = {shift2_joules:.4e} J")
    print(f"   In quecto-electronvolts, the energy shift is: {shift2_qeV:.4e} qeV")

    # Final answer in the required format
    print("\n<<<")
    print(f"The energy shift for the J=0 to J=1 transition is {shift1_qeV:.4e} qeV.")
    print(f"The energy shift for the J=1 to J=2 transition is {shift2_qeV:.4e} qeV.")
    print(">>>")

if __name__ == '__main__':
    calculate_hbr_distortion_shifts()