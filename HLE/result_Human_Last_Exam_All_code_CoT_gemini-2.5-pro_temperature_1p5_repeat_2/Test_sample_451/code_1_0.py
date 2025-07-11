import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.

    This script follows these steps:
    1.  Define physical constants and input parameters for H-Br in SI units.
    2.  Calculate the reduced mass (mu) and moment of inertia (I) of the molecule.
    3.  Calculate the vibrational angular frequency (omega).
    4.  Compute the rotational constant (B) and the centrifugal distortion constant (D) in Joules.
    5.  The energy shift for a transition J -> J+1 is given by 4*D*(J+1)^3.
        This is calculated for the J=0->1 and J=1->2 transitions.
    6.  Convert the final energy shifts from Joules to quecto-electronvolts (qeV) and print the results.
    """
    # 1. Define constants and inputs
    # Given parameters
    m_H_amu = 1.008                   # Mass of Hydrogen in amu
    m_Br_amu = 79.904                 # Mass of Bromine in amu
    r0_pm = 141.4                     # Bond length in pm
    k_Nm = 400                        # Force constant in N/m

    # Physical constants in SI units
    amu_to_kg = 1.66053906660e-27     # kg/amu
    hbar = 1.054571817e-34            # J*s (Reduced Planck constant)
    e_charge = 1.602176634e-19         # C (Elementary charge for eV conversion)

    # Convert inputs to SI units
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    r0_m = r0_pm * 1e-12

    # 2. Calculate molecular properties
    # Reduced mass (mu)
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Moment of inertia (I)
    I = mu * r0_m**2

    # 3. Calculate vibrational frequency
    # Vibrational angular frequency (omega)
    omega = math.sqrt(k_Nm / mu)

    # 4. Compute spectroscopic constants
    # Rotational constant (B) in Joules
    B = hbar**2 / (2 * I)

    # Centrifugal distortion constant (D) in Joules
    D = (4 * B**3) / ((hbar * omega)**2)

    # 5. Calculate energy shifts for each transition
    # For J=0 -> J=1 transition, the initial state is J=0
    shift_J0_to_J1_joules = 4 * D * (0 + 1)**3
    
    # For J=1 -> J=2 transition, the initial state is J=1
    shift_J1_to_J2_joules = 4 * D * (1 + 1)**3

    # 6. Convert to qeV and print
    # 1 qeV = 1e-30 eV = 1e-30 * e_charge Joules
    joules_per_qeV = 1e-30 * e_charge

    shift_J0_to_J1_qeV = shift_J0_to_J1_joules / joules_per_qeV
    shift_J1_to_J2_qeV = shift_J1_to_J2_joules / joules_per_qeV

    print("--- Intermediate Values ---")
    print(f"Reduced mass (mu): {mu:.5e} kg")
    print(f"Moment of inertia (I): {I:.5e} kg*m^2")
    print(f"Centrifugal distortion constant (D): {D:.5e} J")
    print("-" * 30)

    # Output for J=0 -> J=1 transition
    print("Energy shift for rotational transition J=0 -> J=1:")
    print(f"Equation: Shift = 4 * {D:.5e} J * ({0} + 1)^3")
    print(f"Calculated Shift = {shift_J0_to_J1_joules:.5e} J")
    print(f"Shift in qeV = {shift_J0_to_J1_joules:.5e} J / {joules_per_qeV:.5e} J/qeV")
    print(f"The energy shift for the J = 0 to J = 1 transition is {shift_J0_to_J1_qeV:.4e} qeV.\n")

    # Output for J=1 -> J=2 transition
    print("Energy shift for rotational transition J=1 -> J=2:")
    print(f"Equation: Shift = 4 * {D:.5e} J * ({1} + 1)^3")
    print(f"Calculated Shift = {shift_J1_to_J2_joules:.5e} J")
    print(f"Shift in qeV = {shift_J1_to_J2_joules:.5e} J / {joules_per_qeV:.5e} J/qeV")
    print(f"The energy shift for the J = 1 to J = 2 transition is {shift_J1_to_J2_qeV:.4e} qeV.")

calculate_centrifugal_distortion_shifts()