import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br.
    The code follows the plan outlined above, calculating molecular constants
    step-by-step and then determining the energy shifts for the specified
    rotational transitions.
    """

    # --- Configuration: Given and Physical Constants ---

    # Given parameters for H-Br
    r0_pm = 141.4      # Bond length in pm
    k_force = 400.0    # Force constant in N/m
    mH_amu = 1.008     # Mass of Hydrogen in amu
    mBr_amu = 79.904   # Mass of Bromine in amu

    # Physical constants in SI units
    h = 6.62607015e-34      # Planck's constant in J·s
    e_charge = 1.602176634e-19 # Elementary charge in C (for J to eV conversion)
    amu_to_kg = 1.660539e-27   # a.m.u to kg conversion factor

    print("This script calculates the energy shifts due to centrifugal distortion for H-Br.\n")
    print("--- Intermediate Calculations ---\n")

    # Step 1: Calculate Reduced Mass (μ) in kg
    r0_m = r0_pm * 1e-12
    mH_kg = mH_amu * amu_to_kg
    mBr_kg = mBr_amu * amu_to_kg
    mu = (mH_kg * mBr_kg) / (mH_kg + mBr_kg)
    print(f"1. Reduced mass μ = (m_H * m_Br) / (m_H + m_Br)")
    print(f"   μ = ({mH_kg:.4e} kg * {mBr_kg:.4e} kg) / ({mH_kg:.4e} kg + {mBr_kg:.4e} kg) = {mu:.5e} kg\n")

    # Step 2: Calculate Moment of Inertia (I)
    I = mu * r0_m**2
    print(f"2. Moment of inertia I = μ * r₀²")
    print(f"   I = {mu:.5e} kg * ({r0_m:.5e} m)² = {I:.5e} kg·m²\n")

    # Step 3: Calculate Rotational Constant (B) in Joules
    B_energy = h**2 / (8 * math.pi**2 * I)
    print(f"3. Rotational constant B = h² / (8π²I)")
    print(f"   B = ({h:.5e} J·s)² / (8π² * {I:.5e} kg·m²) = {B_energy:.5e} J\n")

    # Step 4: Calculate Fundamental Vibrational Energy (hν) in Joules
    hbar = h / (2 * math.pi)
    omega = math.sqrt(k_force / mu)
    hv_energy = hbar * omega
    print(f"4. Vibrational energy hν = ħ√(k/μ)")
    print(f"   hν = {hbar:.5e} J·s * √({k_force} N/m / {mu:.5e} kg) = {hv_energy:.5e} J\n")

    # Step 5: Calculate Centrifugal Distortion Constant (D) in Joules
    D_energy = 4 * B_energy**3 / hv_energy**2
    print(f"5. Centrifugal distortion constant D = 4B³ / (hν)²")
    print(f"   D = 4 * ({B_energy:.5e} J)³ / ({hv_energy:.5e} J)² = {D_energy:.5e} J\n")
    
    print("-" * 50)
    print("\n--- Final Energy Shift Calculations ---\n")

    # --- Calculation for Transition J=0 -> J=1 ---
    print("For transition J=0 -> J=1:")
    # The energy shift is ΔE_dist = E_dist(J=1) - E_dist(J=0)
    # E_dist(J) = -D*J²*(J+1)². So, E_dist(1) = -D*1²*2² = -4D, and E_dist(0) = 0.
    # ΔE_dist = -4D
    delta_E_J_1 = -4 * D_energy
    delta_E_eV_1 = delta_E_J_1 / e_charge
    delta_E_qeV_1 = delta_E_eV_1 * 1e30
    
    print(f"   Shift equation: ΔE = E_dist(1) - E_dist(0) = -4 * D")
    print(f"   ΔE = -4 * {D_energy:.5e} J = {delta_E_J_1:.5e} J")
    print(f"   Converting to qeV: ({delta_E_J_1:.5e} J / {e_charge:.5e} C) * 10^30 = {delta_E_qeV_1:.5e} qeV\n")

    # --- Calculation for Transition J=1 -> J=2 ---
    print("For transition J=1 -> J=2:")
    # The energy shift is ΔE_dist = E_dist(J=2) - E_dist(J=1)
    # E_dist(2) = -D*2²*3² = -36D, E_dist(1) = -D*1²*2² = -4D.
    # ΔE_dist = -36D - (-4D) = -32D
    delta_E_J_2 = -32 * D_energy
    delta_E_eV_2 = delta_E_J_2 / e_charge
    delta_E_qeV_2 = delta_E_eV_2 * 1e30

    print(f"   Shift equation: ΔE = E_dist(2) - E_dist(1) = -36*D - (-4*D) = -32 * D")
    print(f"   ΔE = -32 * {D_energy:.5e} J = {delta_E_J_2:.5e} J")
    print(f"   Converting to qeV: ({delta_E_J_2:.5e} J / {e_charge:.5e} C) * 10^30 = {delta_E_qeV_2:.5e} qeV\n")

    print("-" * 50)
    print("\n>>> Final Answer Summary <<<")
    print(f"1. Energy shift for J = 0 to J = 1 transition: {delta_E_qeV_1:.4e} qeV")
    print(f"2. Energy shift for J = 1 to J = 2 transition: {delta_E_qeV_2:.4e} qeV")


if __name__ == '__main__':
    calculate_centrifugal_distortion_shifts()