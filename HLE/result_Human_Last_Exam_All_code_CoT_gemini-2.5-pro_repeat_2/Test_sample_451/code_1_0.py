import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- Given Constants ---
    m_H_amu = 1.008         # Atomic mass of Hydrogen (amu)
    m_Br_amu = 79.904       # Atomic mass of Bromine (amu)
    r0_pm = 141.4           # Bond length (pm)
    k_Nm = 400.0            # Force constant (N/m)

    # --- Physical Constants ---
    amu_to_kg = 1.660539e-27
    pm_to_m = 1e-12
    h_bar = 1.054571817e-34  # Reduced Planck constant (J*s)
    e_charge = 1.602176634e-19 # Elementary charge (C)

    # --- Step 1: Calculate reduced mass in kg ---
    mu_amu = (m_H_amu * m_Br_amu) / (m_H_amu + m_Br_amu)
    mu_kg = mu_amu * amu_to_kg

    # --- Step 2: Convert bond length to meters ---
    r0_m = r0_pm * pm_to_m

    # --- Step 3: Calculate the centrifugal distortion constant D in Joules ---
    # Formula: D = ħ⁴ / (2 * k * μ² * r₀⁶)
    D_joules = (h_bar**4) / (2 * k_Nm * mu_kg**2 * r0_m**6)

    # --- Step 4: Define conversion factor from Joules to qeV ---
    # 1 qeV = 1e-30 eV; 1 eV = e_charge J
    J_to_qeV = 1e30 / e_charge

    # --- Step 5: Calculate and print results for each transition ---
    
    # Transition 1: J=0 -> J=1
    J1 = 0
    shift_factor_1 = -4 * (J1 + 1)**3
    E_shift_1_J = shift_factor_1 * D_joules
    E_shift_1_qeV = E_shift_1_J * J_to_qeV

    print("1. Energy shift for J=0 -> J=1 transition:")
    # The final equation with all numbers: ΔE = shift_factor * D = result
    print(f"ΔE = {shift_factor_1} * {D_joules:.4e} J = {E_shift_1_qeV:.4e} qeV")
    print("-" * 40)

    # Transition 2: J=1 -> J=2
    J2 = 1
    shift_factor_2 = -4 * (J2 + 1)**3
    E_shift_2_J = shift_factor_2 * D_joules
    E_shift_2_qeV = E_shift_2_J * J_to_qeV
    
    print("2. Energy shift for J=1 -> J=2 transition:")
    # The final equation with all numbers: ΔE = shift_factor * D = result
    print(f"ΔE = {shift_factor_2} * {D_joules:.4e} J = {E_shift_2_qeV:.4e} qeV")

calculate_centrifugal_distortion_shifts()