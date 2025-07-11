import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for rotational
    transitions in an H-Br molecule.
    """
    # --- Physical Constants ---
    h = 6.62607015e-34      # Planck's constant (J·s)
    amu_to_kg = 1.66053906660e-27 # kg/amu conversion factor
    pm_to_m = 1e-12             # pm/m conversion factor
    J_to_eV = 1.602176634e-19   # J/eV conversion factor
    eV_to_qeV = 1e30            # qeV/eV conversion factor

    # --- Given Values for H-Br ---
    m_H_amu = 1.008             # Mass of Hydrogen (amu)
    m_Br_amu = 79.904           # Mass of Bromine (amu)
    r0_pm = 141.4               # Bond length (pm)
    k_Nm = 400.0                # Force constant (N/m)

    # --- Step 1: Convert given values to SI units ---
    m_H = m_H_amu * amu_to_kg
    m_Br = m_Br_amu * amu_to_kg
    r0 = r0_pm * pm_to_m

    # --- Step 2: Calculate molecular properties ---
    # Reduced mass (μ)
    mu = (m_H * m_Br) / (m_H + m_Br)
    
    # Moment of inertia (I)
    I = mu * r0**2
    
    # Rotational constant (B) in Joules
    B = h**2 / (8 * math.pi**2 * I)
    
    # Vibrational frequency (f_vib) in Hz
    f_vib = (1 / (2 * math.pi)) * math.sqrt(k_Nm / mu)
    
    # Centrifugal distortion constant (D) in Joules
    D = (4 * B**3) / (h**2 * f_vib**2)

    print("This script calculates the energy shifts due to centrifugal distortion for H-Br.")
    print("\n--- 1. Calculated Molecular Properties ---")
    print(f"Reduced mass (μ): {mu:.4e} kg")
    print(f"Moment of inertia (I): {I:.4e} kg·m²")
    print(f"Vibrational frequency (f_vib): {f_vib:.4e} Hz")
    print(f"Rotational constant (B): {B:.4e} J")
    print(f"Centrifugal distortion constant (D): {D:.4e} J")

    # --- Step 3: Calculate Energy Shift for J=0 -> J=1 ---
    print("\n--- 2. Energy Shift for J=0 to J=1 transition ---")
    J_initial_1 = 0
    print("The energy shift (ΔE_dist) for a transition from J to J+1 is given by the formula: ΔE_dist = -4 * D * (J + 1)³")
    
    # Calculation
    delta_E_dist_J1 = -4 * D * (J_initial_1 + 1)**3
    delta_E_dist_qeV1 = (delta_E_dist_J1 / J_to_eV) * eV_to_qeV

    print(f"For the J={J_initial_1} to J={J_initial_1 + 1} transition, the equation is:")
    print(f"ΔE_dist = -4 * ({D:.4e} J) * ({J_initial_1} + 1)³ = {delta_E_dist_J1:.4e} J")
    print(f"Energy shift in qeV = ({delta_E_dist_J1:.4e} J / {J_to_eV:.4e} J/eV) * {eV_to_qeV:.0e} qeV/eV")
    print(f"Final Energy Shift for J=0 -> J=1: {delta_E_dist_qeV1:.4e} qeV")

    # --- Step 4: Calculate Energy Shift for J=1 -> J=2 ---
    print("\n--- 3. Energy Shift for J=1 to J=2 transition ---")
    J_initial_2 = 1
    print("The energy shift (ΔE_dist) for a transition from J to J+1 is given by the formula: ΔE_dist = -4 * D * (J + 1)³")

    # Calculation
    delta_E_dist_J2 = -4 * D * (J_initial_2 + 1)**3
    delta_E_dist_qeV2 = (delta_E_dist_J2 / J_to_eV) * eV_to_qeV

    print(f"For the J={J_initial_2} to J={J_initial_2 + 1} transition, the equation is:")
    print(f"ΔE_dist = -4 * ({D:.4e} J) * ({J_initial_2} + 1)³ = {delta_E_dist_J2:.4e} J")
    print(f"Energy shift in qeV = ({delta_E_dist_J2:.4e} J / {J_to_eV:.4e} J/eV) * {eV_to_qeV:.0e} qeV/eV")
    print(f"Final Energy Shift for J=1 -> J=2: {delta_E_dist_qeV2:.4e} qeV")
    
    # Returning final answer for the wrapper
    return f"{delta_E_dist_qeV1:.4e}, {delta_E_dist_qeV2:.4e}"

# Execute the calculation and print the final answer in the required format
final_answers = calculate_centrifugal_distortion_shifts()
print(f"\n<<<{final_answers}>>>")
