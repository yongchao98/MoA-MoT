import math

def calculate_centrifugal_distortion_shift():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- Step 1: Define Physical Constants and Input Parameters ---
    h_bar = 1.054571817e-34      # Reduced Planck constant (J·s)
    amu_to_kg = 1.66053906660e-27  # Atomic mass unit to kg
    e_charge = 1.602176634e-19     # Elementary charge (C) for J to eV conversion
    pm_to_m = 1e-12                # Picometer to meter conversion

    # Given data for H-Br
    r0_pm = 141.4  # Bond length in pm
    k = 400.0      # Force constant in N/m
    m_H_amu = 1.008  # Mass of Hydrogen in amu
    m_Br_amu = 79.904 # Mass of Bromine in amu

    # --- Step 2: Calculate Molecular Properties in SI units ---
    r0 = r0_pm * pm_to_m
    m_H = m_H_amu * amu_to_kg
    m_Br = m_Br_amu * amu_to_kg

    # Calculate reduced mass (mu)
    mu = (m_H * m_Br) / (m_H + m_Br)

    # --- Step 3: Calculate the Centrifugal Distortion Constant (D) in Joules ---
    # The formula is D = ħ⁴ / (2 * k * μ² * r₀⁶)
    D_J = (h_bar**4) / (2 * k * mu**2 * r0**6)

    # --- Step 4 & 5: Calculate Energy Shift for Each Transition ---
    # The shift in transition energy from J -> J+1 is ΔE = -4 * D * (J_initial + 1)³
    
    # Transition 1: J = 0 -> J = 1
    J1_initial = 0
    delta_E1_J = -4 * D_J * (J1_initial + 1)**3

    # Transition 2: J = 1 -> J = 2
    J2_initial = 1
    delta_E2_J = -4 * D_J * (J2_initial + 1)**3

    # --- Step 6: Convert Energies to quecto-electronvolts (qeV) ---
    # Conversion factor from Joules to qeV
    joules_to_qeV = 1e30 / e_charge
    
    delta_E1_qeV = delta_E1_J * joules_to_qeV
    delta_E2_qeV = delta_E2_J * joules_to_qeV

    # --- Final Output ---
    print("--- Calculation for J=0 -> J=1 transition ---")
    print(f"The energy shift due to centrifugal distortion is given by: ΔE = -4 * D * (J_initial + 1)³")
    print(f"Centrifugal distortion constant, D = {D_J:.3e} J")
    print(f"ΔE = -4 * ({D_J:.3e} J) * ({J1_initial} + 1)³ = {delta_E1_J:.3e} J")
    print(f"Converting to quecto-electronvolts (qeV): {delta_E1_J:.3e} J * ({joules_to_qeV:.3e} qeV/J)")
    print(f"Energy shift for J=0 -> J=1: {delta_E1_qeV:.3e} qeV\n")

    print("--- Calculation for J=1 -> J=2 transition ---")
    print(f"The energy shift due to centrifugal distortion is given by: ΔE = -4 * D * (J_initial + 1)³")
    print(f"Centrifugal distortion constant, D = {D_J:.3e} J")
    print(f"ΔE = -4 * ({D_J:.3e} J) * ({J2_initial} + 1)³ = {delta_E2_J:.3e} J")
    print(f"Converting to quecto-electronvolts (qeV): {delta_E2_J:.3e} J * ({joules_to_qeV:.3e} qeV/J)")
    print(f"Energy shift for J=1 -> J=2: {delta_E2_qeV:.3e} qeV")

if __name__ == '__main__':
    calculate_centrifugal_distortion_shift()