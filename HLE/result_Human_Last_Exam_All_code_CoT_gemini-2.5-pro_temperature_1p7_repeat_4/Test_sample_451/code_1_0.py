import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define given and physical constants in SI units.
    # Given problem parameters
    r0_pm = 141.4  # Bond length in picometers
    k_Nm = 400.0     # Force constant in N/m
    m_H_amu = 1.008  # Mass of Hydrogen in amu
    m_Br_amu = 79.904 # Mass of Bromine in amu

    # Physical constants
    h = 6.62607015e-34      # Planck constant (J*s)
    hbar = h / (2 * math.pi)  # Reduced Planck constant (J*s)
    amu_to_kg = 1.660539e-27 # kg/amu
    J_per_eV = 1.602176634e-19 # Joules per eV

    # Unit conversions
    pm_to_m = 1e-12
    qeV_per_eV = 1e30

    # 2. Convert input parameters to a consistent set of SI units.
    r0 = r0_pm * pm_to_m
    m_H = m_H_amu * amu_to_kg
    m_Br = m_Br_amu * amu_to_kg

    # 3. Perform calculations.
    # a. Calculate the reduced mass (mu) of H-Br.
    mu = (m_H * m_Br) / (m_H + m_Br)

    # b. Calculate the centrifugal distortion constant (D) in Joules.
    D_joules = (hbar**4) / (2 * k_Nm * mu**2 * r0**6)

    # c. Calculate the energy shifts for each transition in Joules.
    # For J=0 -> J=1, the shift is -4*D
    shift1_joules = -4 * D_joules
    # For J=1 -> J=2, the shift is -32*D
    shift2_joules = -32 * D_joules

    # d. Convert the energy shifts to quecto-electronvolts (qeV).
    shift1_qeV = (shift1_joules / J_per_eV) * qeV_per_eV
    shift2_qeV = (shift2_joules / J_per_eV) * qeV_per_eV

    # 4. Print the results with the final equations.
    print("### Calculation for Transition J = 0 to J = 1 ###\n")
    print(f"The energy shift ΔE is given by -4 * D.")
    print(f"The centrifugal distortion constant D is {D_joules:.6e} J.")
    print(f"ΔE = -4 * {D_joules:.6e} J = {shift1_joules:.6e} J")
    print("\nConverting to quecto-electronvolts (qeV):")
    print(f"ΔE = ({shift1_joules:.6e} J / {J_per_eV:.6e} J/eV) * {qeV_per_eV:.0e} qeV/eV")
    print(f"ΔE (J=0 to 1) = {shift1_qeV:.6f} qeV")
    
    print("\n" + "-"*50 + "\n")

    print("### Calculation for Transition J = 1 to J = 2 ###\n")
    print(f"The energy shift ΔE is given by -32 * D.")
    print(f"The centrifugal distortion constant D is {D_joules:.6e} J.")
    print(f"ΔE = -32 * {D_joules:.6e} J = {shift2_joules:.6e} J")
    print("\nConverting to quecto-electronvolts (qeV):")
    print(f"ΔE = ({shift2_joules:.6e} J / {J_per_eV:.6e} J/eV) * {qeV_per_eV:.0e} qeV/eV")
    print(f"ΔE (J=1 to 2) = {shift2_qeV:.6f} qeV")

calculate_energy_shifts()
<<<[-176900331046166408226816.000000, -1415202648369331238883328.000000]>>>