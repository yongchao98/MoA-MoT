import scipy.constants as const

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define Constants
    r0 = 141.4e-12  # bond length in meters
    k = 400.0  # force constant in N/m
    m_H_amu = 1.008  # mass of Hydrogen in amu
    m_Br_amu = 79.904 # mass of Bromine in amu
    
    # Physical constants from scipy
    amu_to_kg = const.physical_constants['atomic mass constant'][0]
    hbar = const.hbar
    e = const.e

    # Unit conversions
    # 1 qeV = 1e-30 eV
    # 1 eV = e Joules
    joule_to_qeV = 1 / (e * 1e-30)

    # 2. Calculate Reduced Mass (mu) in kg
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # 3. Calculate Centrifugal Distortion Constant (D_J) in Joules
    # Formula: D_J = hbar^4 / (2 * k * mu^2 * r0^6)
    D_J = (hbar**4) / (2 * k * mu**2 * r0**6)

    # 4. Calculate Energy Shifts for Transitions
    # Formula: dE_shift = -4 * D_J * (J+1)^3

    # Transition 1: J = 0 -> J = 1
    J1 = 0
    shift1_J = -4 * D_J * (J1 + 1)**3
    shift1_qeV = shift1_J * joule_to_qeV
    
    # Transition 2: J = 1 -> J = 2
    J2 = 1
    shift2_J = -4 * D_J * (J2 + 1)**3
    shift2_qeV = shift2_J * joule_to_qeV
    
    # Print the results as requested
    print(f"The calculated centrifugal distortion constant D_J = {D_J:.4e} J.\n")

    print("--- Transition 1: J = 0 to J = 1 ---")
    print(f"Energy Shift Equation: ΔE = -4 * D_J * (J+1)^3")
    print(f"ΔE = -4 * ({D_J:.4e} J) * ({J1}+1)^3")
    print(f"ΔE = {shift1_J:.4e} J")
    print(f"Energy shift in quecto-electronvolts (qeV): {shift1_qeV:.4e} qeV\n")

    print("--- Transition 2: J = 1 to J = 2 ---")
    print(f"Energy Shift Equation: ΔE = -4 * D_J * (J+1)^3")
    print(f"ΔE = -4 * ({D_J:.4e} J) * ({J2}+1)^3")
    print(f"ΔE = {shift2_J:.4e} J")
    print(f"Energy shift in quecto-electronvolts (qeV): {shift2_qeV:.4e} qeV")
    
    # The final answer block is expected to be parsed programmatically.
    # It will contain the two requested numerical answers.
    global final_answer
    final_answer = f"<<<{shift1_qeV:.4e}, {shift2_qeV:.4e}>>>"

# Execute the calculation
calculate_energy_shifts()
# This print statement is for the final answer block as per instructions.
print(final_answer)