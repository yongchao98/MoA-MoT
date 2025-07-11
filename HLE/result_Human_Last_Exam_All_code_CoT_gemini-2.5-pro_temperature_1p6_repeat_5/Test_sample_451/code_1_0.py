import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br.
    """
    # 1. Define Constants
    h = 6.62607015e-34  # Planck's constant (J·s)
    hbar = h / (2 * math.pi)  # Reduced Planck's constant (J·s)
    amu_to_kg = 1.66053906660e-27  # Atomic mass unit to kg conversion
    pm_to_m = 1e-12  # Picometer to meter conversion
    eV_to_J = 1.602176634e-19  # Electron-volt to Joules conversion
    qeV_to_J = 1e-30 * eV_to_J # quecto-eV to Joules conversion

    # Input parameters for H-Br
    r0_pm = 141.4  # Bond length (pm)
    k = 400.0  # Force constant (N/m)
    m_H_amu = 1.008  # Mass of Hydrogen (amu)
    m_Br_amu = 79.904 # Mass of Bromine (amu)

    # Convert inputs to SI units
    r0 = r0_pm * pm_to_m
    m_H = m_H_amu * amu_to_kg
    m_Br = m_Br_amu * amu_to_kg

    # 2. Calculate Molecular Properties
    # Reduced mass (mu)
    mu = (m_H * m_Br) / (m_H + m_Br)

    # Moment of inertia (I)
    I = mu * r0**2

    # 3. Calculate Spectroscopic Constants
    # Rotational constant (B) in Joules
    B_joule = hbar**2 / (2 * I)

    # Angular vibrational frequency (omega) in rad/s
    omega_rad_s = math.sqrt(k / mu)

    # Energy of a vibrational quantum (hw) in Joules
    hw_joule = hbar * omega_rad_s
    
    # Centrifugal distortion constant (D) in Joules
    D_joule = (4 * B_joule**3) / (hw_joule**2)

    # 4. Calculate Energy Shifts in Joules
    # Transition J=0 -> J=1 (lower state J=0)
    J1 = 0
    shift1_joule = -4 * D_joule * (J1 + 1)**3
    
    # Transition J=1 -> J=2 (lower state J=1)
    J2 = 1
    shift2_joule = -4 * D_joule * (J2 + 1)**3

    # 5. Convert and Output Results
    # Convert from Joules to qeV
    shift1_qeV = shift1_joule / qeV_to_J
    shift2_qeV = shift2_joule / qeV_to_J

    print(f"The calculated centrifugal distortion constant is D = {D_joule:.4e} J.")
    print("-" * 50)
    
    print("1. For the rotational transition from J = 0 to J = 1:")
    # The final equation for the shift is: ΔE = -4 * D * (J+1)³
    print(f"   Equation: ΔE = -4 * ({D_joule:.4e} J) * ({J1}+1)³")
    print(f"   The energy shift is {shift1_joule:.4e} J.")
    print(f"   In quecto-electronvolts, ΔE = {shift1_qeV:.4e} qeV.")
    
    print("\n" + "-" * 50)

    print("2. For the rotational transition from J = 1 to J = 2:")
    # The final equation for the shift is: ΔE = -4 * D * (J+1)³
    print(f"   Equation: ΔE = -4 * ({D_joule:.4e} J) * ({J2}+1)³")
    print(f"   The energy shift is {shift2_joule:.4e} J.")
    print(f"   In quecto-electronvolts, ΔE = {shift2_qeV:.4e} qeV.")

if __name__ == '__main__':
    calculate_energy_shifts()