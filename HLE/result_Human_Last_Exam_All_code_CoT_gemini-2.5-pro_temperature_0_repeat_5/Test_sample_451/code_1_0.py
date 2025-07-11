import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br
    for J=0->1 and J=1->2 transitions.
    """
    # --- Constants ---
    h = 6.62607015e-34      # Planck's constant in J·s
    hbar = h / (2 * math.pi)  # Reduced Planck's constant in J·s
    amu_to_kg = 1.660539e-27  # AMU to kg conversion factor
    e_charge = 1.602176634e-19 # Elementary charge in Coulombs

    # --- Given values for H-Br ---
    r0 = 141.4e-12          # Bond length in meters
    k = 400                 # Force constant in N/m
    m_H_amu = 1.008         # Mass of Hydrogen in amu
    m_Br_amu = 79.904       # Mass of Bromine in amu

    # --- Step 1: Calculate reduced mass (mu) in kg ---
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # --- Step 2: Calculate moment of inertia (I) ---
    I = mu * r0**2

    # --- Step 3: Calculate rotational constant (B) and vibrational frequency (omega) ---
    # Rotational constant B in Joules
    B_joules = hbar**2 / (2 * I)
    # Angular vibrational frequency omega in rad/s
    omega = math.sqrt(k / mu)

    # --- Step 4: Calculate centrifugal distortion constant (D) in Joules ---
    # The formula is D = 4 * B^3 / (hbar^2 * omega^2)
    D_joules = 4 * B_joules**3 / (hbar**2 * omega**2)

    # --- Step 5: Calculate energy shifts for transitions in Joules ---
    # For J=0 -> J=1, the shift is Δ(ΔE) = E_dist(J=1) - E_dist(J=0) = -D*[1^2*(2^2)] - 0 = -4*D
    shift_0_to_1_joules = -4 * D_joules
    
    # For J=1 -> J=2, the shift is Δ(ΔE) = E_dist(J=2) - E_dist(J=1) = -D*[2^2*(3^2)] - [-D*1^2*(2^2)] = -36*D + 4*D = -32*D
    shift_1_to_2_joules = -32 * D_joules

    # --- Step 6: Convert energy shifts to quecto-electronvolts (qeV) ---
    # 1 qeV = 1e-30 eV
    # 1 eV = 1.602176634e-19 J
    # So, 1 qeV = 1e-30 * 1.602176634e-19 J = 1.602176634e-49 J
    joules_to_qeV = 1 / (e_charge * 1e-30)

    shift_0_to_1_qeV = shift_0_to_1_joules * joules_to_qeV
    shift_1_to_2_qeV = shift_1_to_2_joules * joules_to_qeV

    # --- Step 7: Print the results ---
    print("Calculation of Centrifugal Distortion Energy Shifts for H-Br\n")
    print(f"The calculated centrifugal distortion constant is D = {D_joules:.4e} J.\n")
    
    print("1. For the rotational transition from J = 0 to J = 1:")
    print(f"   The energy shift is given by Δ(ΔE) = -4 * D")
    print(f"   Δ(ΔE) = -4 * {D_joules:.4e} J = {shift_0_to_1_joules:.4e} J")
    print(f"   In quecto-electronvolts, this is:")
    print(f"   Δ(ΔE) = {shift_0_to_1_qeV:.4f} qeV\n")

    print("2. For the rotational transition from J = 1 to J = 2:")
    print(f"   The energy shift is given by Δ(ΔE) = -32 * D")
    print(f"   Δ(ΔE) = -32 * {D_joules:.4e} J = {shift_1_to_2_joules:.4e} J")
    print(f"   In quecto-electronvolts, this is:")
    print(f"   Δ(ΔE) = {shift_1_to_2_qeV:.4f} qeV")

if __name__ == '__main__':
    calculate_centrifugal_distortion_shifts()