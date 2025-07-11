import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define Constants
    # Given values
    r0_pm = 141.4  # bond length in picometers
    k_Nm = 400.0   # force constant in N/m
    m_H_amu = 1.008  # mass of Hydrogen in amu
    m_Br_amu = 79.904 # mass of Bromine in amu

    # Fundamental physical constants
    amu_to_kg = 1.660539e-27 # conversion from amu to kg
    pm_to_m = 1e-12          # conversion from pm to m
    hbar = 1.054571817e-34   # Reduced Planck's constant in J*s
    e_charge = 1.602176634e-19 # Elementary charge in Coulombs

    # Unit conversions
    r0_m = r0_pm * pm_to_m
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg

    # 2. Calculate Molecular Properties
    # Calculate reduced mass (mu) in kg
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Calculate moment of inertia (I) in kg*m^2
    I = mu_kg * r0_m**2

    # 3. Calculate Spectroscopic Constants
    # Calculate rotational constant (B) in Joules
    B_joules = hbar**2 / (2 * I)

    # Calculate angular vibrational frequency (omega_e) in rad/s
    omega_e = math.sqrt(k_Nm / mu_kg)
    # Calculate vibrational energy quantum (E_vib) in Joules
    E_vib_joules = hbar * omega_e

    # Calculate centrifugal distortion constant (D) in Joules
    D_joules = 4 * (B_joules**3) / (E_vib_joules**2)

    # 4. Calculate Energy Shifts for each transition
    transitions = [0, 1]  # J values for the initial states

    print(f"Calculated centrifugal distortion constant D = {D_joules:.4e} J\n")

    for J in transitions:
        # Calculate energy shift in Joules
        delta_E_J = -4 * D_joules * (J + 1)**3

        # 5. Convert Units
        # Convert energy shift to electronvolts (eV)
        delta_E_eV = delta_E_J / e_charge
        # Convert energy shift to quecto-electronvolts (qeV)
        delta_E_qeV = delta_E_eV * 1e30

        # Print the final results for the current transition
        print(f"For the transition from J = {J} to J = {J+1}:")
        print(f"The energy shift (ΔE) is calculated as: ΔE = -4 * D * (J+1)³")
        print(f"ΔE = -4 * {D_joules:.4e} J * ({J}+1)³ = {delta_E_J:.4e} J")
        print(f"In quecto-electronvolts, this is: {delta_E_qeV:.4e} qeV\n")

# Run the calculation and print the results
calculate_energy_shifts()

# <<<
# For the transition from J = 0 to J = 1:
# The energy shift (ΔE) is calculated as: ΔE = -4 * D * (J+1)³
# ΔE = -4 * 7.0787e-27 J * (0+1)³ = -2.8315e-26 J
# In quecto-electronvolts, this is: -1.7673e+23 qeV
#
# For the transition from J = 1 to J = 2:
# The energy shift (ΔE) is calculated as: ΔE = -4 * D * (J+1)³
# ΔE = -4 * 7.0787e-27 J * (1+1)³ = -2.2652e-25 J
# In quecto-electronvolts, this is: -1.4138e+24 qeV
# >>>