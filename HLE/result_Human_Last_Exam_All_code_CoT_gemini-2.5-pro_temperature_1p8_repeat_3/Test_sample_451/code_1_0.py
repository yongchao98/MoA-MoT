import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define Physical Constants and Input Parameters (SI units)
    H = 6.62607015e-34  # Planck's constant (J*s)
    HBAR = H / (2 * math.pi) # Reduced Planck's constant (J*s)
    AMU = 1.66053906660e-27 # Atomic mass unit (kg)
    E_CHARGE = 1.602176634e-19 # Elementary charge (C)

    # Given parameters for H-Br
    r0 = 141.4e-12      # Bond length (m)
    k_force = 400      # Force constant (N/m)
    m_h = 1.008 * AMU    # Mass of Hydrogen (kg)
    m_br = 79.904 * AMU  # Mass of Bromine (kg)

    # 2. Calculate Molecular Properties
    # Reduced Mass (mu)
    mu = (m_h * m_br) / (m_h + m_br)

    # Moment of Inertia (I)
    I = mu * r0**2

    # Vibrational Energy (E_vib)
    omega = math.sqrt(k_force / mu) # Angular frequency in rad/s
    E_vib = HBAR * omega

    # 3. Calculate Spectroscopic Constants (in Joules)
    # Rotational Constant (B)
    B_joule = HBAR**2 / (2 * I)

    # Centrifugal Distortion Constant (D)
    D_joule = (4 * B_joule**3) / (E_vib**2)

    # 4. Convert D to qeV
    D_eV = D_joule / E_CHARGE
    D_qeV = D_eV * 1e30

    print(f"Calculated centrifugal distortion constant D = {D_qeV:.4f} qeV")
    print("-" * 20)
    
    # 5. Calculate and print the energy shift for the first transition
    J1 = 0
    j_term1 = (J1 + 1)**3
    shift1_qeV = -4 * D_qeV * j_term1
    print(f"Energy shift for transition J = {J1} -> J = {J1+1}:")
    print(f"ΔE = -4 * D * (J+1)^3")
    print(f"ΔE = -4 * {D_qeV:.4f} qeV * ({J1}+1)^3 = {shift1_qeV:.4e} qeV")
    print("-" * 20)

    # 6. Calculate and print the energy shift for the second transition
    J2 = 1
    j_term2 = (J2 + 1)**3
    shift2_qeV = -4 * D_qeV * j_term2
    print(f"Energy shift for transition J = {J2} -> J = {J2+1}:")
    print(f"ΔE = -4 * D * (J+1)^3")
    print(f"ΔE = -4 * {D_qeV:.4f} qeV * ({J2}+1)^3 = {shift2_qeV:.4e} qeV")

if __name__ == '__main__':
    calculate_energy_shifts()