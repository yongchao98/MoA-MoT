import math

def calculate_hbr_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define constants
    m_H_amu = 1.008           # amu
    m_Br_amu = 79.904         # amu
    r0_pm = 141.4             # picometers
    k_Nm = 400.0              # N/m
    
    # Physical constants
    amu_to_kg = 1.66053906660e-27  # kg/amu
    pm_to_m = 1.0e-12             # m/pm
    h = 6.62607015e-34            # Planck constant in J·s
    c_m_s = 299792458             # Speed of light in m/s
    c_cm_s = c_m_s * 100          # Speed of light in cm/s
    e_charge = 1.602176634e-19    # Elementary charge in C (for J to eV conversion)

    # 2. Calculate molecular properties
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)
    r0_m = r0_pm * pm_to_m
    I_kg_m2 = mu_kg * r0_m**2

    # 3. Calculate spectroscopic constants in cm^-1
    # Rotational constant B_tilde in cm^-1
    B_tilde_cm1 = h / (8 * math.pi**2 * c_cm_s * I_kg_m2)
    
    # Vibrational frequency nu_e_tilde in cm^-1
    nu_e_tilde_cm1 = (1 / (2 * math.pi * c_cm_s)) * math.sqrt(k_Nm / mu_kg)
    
    # Centrifugal distortion constant D_tilde in cm^-1
    D_tilde_cm1 = (4 * B_tilde_cm1**3) / (nu_e_tilde_cm1**2)

    # 4. Convert D to Joules for energy calculations
    cm1_to_J = h * c_cm_s
    D_joules = D_tilde_cm1 * cm1_to_J

    print(f"Calculated centrifugal distortion constant D = {D_joules:.4e} J\n")

    # 5. Calculate energy shifts for each transition
    J_to_qeV = (1 / e_charge) * 1e30

    # For J=0 -> J=1 transition
    J0 = 0
    shift_J1_joules = -4 * D_joules * (J0 + 1)**3
    shift_J1_qeV = shift_J1_joules * J_to_qeV
    
    print("1. For the J=0 -> J=1 transition:")
    print(f"   ΔE = -4 * D * (J+1)³")
    print(f"   ΔE = -4 * ({D_joules:.4e} J) * ({J0}+1)³")
    print(f"   Energy Shift = {shift_J1_qeV:.4e} qeV")
    print("-" * 20)

    # For J=1 -> J=2 transition
    J1 = 1
    shift_J2_joules = -4 * D_joules * (J1 + 1)**3
    shift_J2_qeV = shift_J2_joules * J_to_qeV

    print("2. For the J=1 -> J=2 transition:")
    print(f"   ΔE = -4 * D * (J+1)³")
    print(f"   ΔE = -4 * ({D_joules:.4e} J) * ({J1}+1)³")
    print(f"   Energy Shift = {shift_J2_qeV:.4e} qeV")

if __name__ == "__main__":
    calculate_hbr_distortion_shifts()

<<<[-1.7763e+23, -1.4210e+24]>>>