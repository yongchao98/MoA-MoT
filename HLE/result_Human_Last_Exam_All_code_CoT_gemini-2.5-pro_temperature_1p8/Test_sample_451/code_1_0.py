import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define Physical Constants
    r_0 = 141.4e-12  # bond length in meters
    k = 400.0  # force constant in N/m
    m_H_amu = 1.008  # mass of hydrogen in amu
    m_Br_amu = 79.904  # mass of bromine in amu

    # Physical constants in SI units
    AMU_TO_KG = 1.660539e-27
    H = 6.62607015e-34  # Planck's constant in J·s
    HBAR = H / (2 * math.pi)  # Reduced Planck's constant in J·s
    EV_TO_J = 1.602176634e-19 # aV to J conversion factor
    QEV_TO_EV = 1e-30 # qeV to eV conversion factor
    J_TO_QEV = 1 / (EV_TO_J * QEV_TO_EV)

    # 2. Calculate Molecular Properties
    # Convert masses to kg
    m_H_kg = m_H_amu * AMU_TO_KG
    m_Br_kg = m_Br_amu * AMU_TO_KG

    # Calculate reduced mass (mu) in kg
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Calculate rotational constant B in Joules
    I = mu * r_0**2  # Moment of inertia
    B_J = HBAR**2 / (2 * I)

    # Calculate centrifugal distortion constant D in Joules
    omega_sq = k / mu  # Angular frequency squared
    D_J = 4 * (B_J**3) / ((HBAR**2) * omega_sq)

    # 3. Calculate Energy Shifts
    # For J=0 -> J=1 transition
    J_i1 = 0
    J_f1 = 1
    # The shift is -D * [J_f^2(J_f+1)^2 - J_i^2(J_i+1)^2] = -D * [1*4 - 0] = -4*D
    delta_E1_J = -4 * D_J

    # For J=1 -> J=2 transition
    J_i2 = 1
    J_f2 = 2
    # The shift is -D * [J_f^2(J_f+1)^2 - J_i^2(J_i+1)^2] = -D * [4*9 - 1*4] = -D * [36-4] = -32*D
    delta_E2_J = -32 * D_J

    # 4. Convert Units and Print Results
    delta_E1_qeV = delta_E1_J * J_TO_QEV
    delta_E2_qeV = delta_E2_J * J_TO_QEV

    print(f"The calculated centrifugal distortion constant, D = {D_J:.4e} J.")
    print("-" * 50)
    
    print("1. For the rotational transition from J = 0 to J = 1:")
    # Using 'g' for general format to avoid scientific notation if not needed,
    # but given the magnitude, it will be used. '.3g' gives 3 significant figures.
    print(f"   ΔE = -4 * D = -4 * {D_J:.4e} J = {delta_E1_J:.4e} J")
    print(f"   ΔE = {delta_E1_qeV:.3e} qeV")
    
    print("\n2. For the rotational transition from J = 1 to J = 2:")
    print(f"   ΔE = -32 * D = -32 * {D_J:.4e} J = {delta_E2_J:.4e} J")
    print(f"   ΔE = {delta_E2_qeV:.3e} qeV")

calculate_energy_shifts()