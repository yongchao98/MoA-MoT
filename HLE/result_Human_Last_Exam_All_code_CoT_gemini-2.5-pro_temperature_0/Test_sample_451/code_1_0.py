import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define Physical Constants and Input Parameters
    # Physical constants
    h = 6.62607015e-34  # Planck's constant (J·s)
    hbar = h / (2 * math.pi)  # Reduced Planck's constant (J·s)
    amu_to_kg = 1.66053906660e-27  # Atomic mass unit to kg
    e = 1.602176634e-19  # Elementary charge (C)

    # Given parameters for H-Br
    r0_pm = 141.4  # Bond length (pm)
    k = 400.0  # Force constant (N/m)
    m_H_amu = 1.008  # Mass of Hydrogen (amu)
    m_Br_amu = 79.904  # Mass of Bromine (amu)

    # Convert parameters to SI units
    r0 = r0_pm * 1e-12  # Bond length (m)
    m_H = m_H_amu * amu_to_kg  # Mass of H (kg)
    m_Br = m_Br_amu * amu_to_kg  # Mass of Br (kg)

    # 2. Calculate Reduced Mass (μ)
    mu = (m_H * m_Br) / (m_H + m_Br)

    # 3. Calculate the Centrifugal Distortion Constant (D)
    # Moment of inertia (I)
    I = mu * r0**2
    # Rotational constant (B) in Joules
    B_J = hbar**2 / (2 * I)
    # Vibrational angular frequency (ω)
    omega = math.sqrt(k / mu)
    # Centrifugal distortion constant (D) in Joules
    D_J = 4 * B_J**3 / (hbar * omega)**2

    # 4. Calculate the Energy Shift for Transitions
    # For J=0 -> J=1, the shift is -4D
    J0 = 0
    delta_E_J01_J = -4 * D_J * (J0 + 1)**3

    # For J=1 -> J=2, the shift is -4D(1+1)^3 = -32D
    J1 = 1
    delta_E_J12_J = -4 * D_J * (J1 + 1)**3

    # 5. Convert to quecto-electronvolts (qeV)
    J_per_eV = e
    J_per_qeV = 1e-30 * J_per_eV
    delta_E_J01_qeV = delta_E_J01_J / J_per_qeV
    delta_E_J12_qeV = delta_E_J12_J / J_per_qeV

    # 6. Format the Output
    print("Calculation Details:")
    print(f"Reduced mass (μ) = {mu:.4e} kg")
    print(f"Centrifugal distortion constant (D) = {D_J:.4e} J\n")

    print("1. For the transition J = 0 -> J = 1:")
    print(f"The energy shift is ΔE = -4 * D * (0+1)³")
    print(f"ΔE = -4 * {D_J:.4e} J = {delta_E_J01_J:.4e} J")
    print(f"ΔE = {delta_E_J01_qeV:.4e} qeV\n")

    print("2. For the transition J = 1 -> J = 2:")
    print(f"The energy shift is ΔE = -4 * D * (1+1)³ = -32 * D")
    print(f"ΔE = -32 * {D_J:.4e} J = {delta_E_J12_J:.4e} J")
    print(f"ΔE = {delta_E_J12_qeV:.4e} qeV")
    
    # Store final answers for the specified format
    global final_answer
    final_answer = f"{delta_E_J01_qeV:.3e}, {delta_E_J12_qeV:.3e}"


# Run the calculation
calculate_energy_shifts()
print(f"\n<<<{final_answer}>>>")