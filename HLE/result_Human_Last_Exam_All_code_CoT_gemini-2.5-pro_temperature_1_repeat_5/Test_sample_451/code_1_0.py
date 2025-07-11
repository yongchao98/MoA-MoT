import scipy.constants as const

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define given parameters and physical constants
    m_H_amu = 1.008      # Atomic mass of Hydrogen (amu)
    m_Br_amu = 79.904    # Atomic mass of Bromine (amu)
    r0_pm = 141.4        # Bond length (pm)
    k_Nm = 400.0         # Force constant (N/m)

    # High-precision physical constants
    amu_to_kg = const.physical_constants['atomic mass unit-kilogram relationship'][0]
    pm_to_m = 1e-12
    hbar = const.hbar
    J_to_eV = 1 / const.e

    # 2. Calculate reduced mass (mu) in kg
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # 3. Calculate moment of inertia (I) in kg*m^2
    r0_m = r0_pm * pm_to_m
    I = mu_kg * r0_m**2

    # 4. Calculate rotational constant (B) in Joules
    B_J = hbar**2 / (2 * I)

    # 5. Calculate the square of the vibrational energy quantum (hbar * omega_e)^2
    # The formula relates it to the force constant k and reduced mass mu
    hbar_omega_e_sq = hbar**2 * k_Nm / mu_kg

    # 6. Calculate centrifugal distortion constant (D) in Joules
    D_J = 4 * B_J**3 / hbar_omega_e_sq

    # 7. Calculate energy shift for the J=0 -> J=1 transition
    # The change in distortion energy for a J -> J+1 transition is -4*D*(J+1)^3
    J0 = 0
    delta_E1_J = -4 * D_J * (J0 + 1)**3

    # 8. Calculate energy shift for the J=1 -> J=2 transition
    J1 = 1
    delta_E2_J = -4 * D_J * (J1 + 1)**3

    # 9. Convert energy shifts from Joules to quecto-electronvolts (qeV)
    # 1 qeV = 1e-30 eV
    J_to_qeV = J_to_eV * 1e30
    delta_E1_qeV = delta_E1_J * J_to_qeV
    delta_E2_qeV = delta_E2_J * J_to_qeV

    # 10. Print the results clearly, showing the final equations
    print("1. Energy shift for rotational transition from J = 0 to J = 1:")
    print(f"ΔE = -4 * D * (0 + 1)^3 = -4 * ({D_J:.4e} J) = {delta_E1_qeV:.4e} qeV")
    print("\n2. Energy shift for rotational transition from J = 1 to J = 2:")
    print(f"ΔE = -4 * D * (1 + 1)^3 = -32 * ({D_J:.4e} J) = {delta_E2_qeV:.4e} qeV")

calculate_energy_shifts()