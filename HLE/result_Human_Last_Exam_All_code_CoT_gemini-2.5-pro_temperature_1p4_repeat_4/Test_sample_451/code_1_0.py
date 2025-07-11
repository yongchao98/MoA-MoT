import math

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- 1. Constants and Given Data ---
    m_H_amu = 1.008  # amu
    m_Br_amu = 79.904 # amu
    r0 = 141.4e-12  # pm to m
    k = 400.0  # N/m
    
    # Physical constants in SI units
    AMU_TO_KG = 1.66053906660e-27
    H_BAR = 1.054571817e-34  # J*s (reduced Planck constant)
    E_CHARGE = 1.602176634e-19 # C (electron charge for J to eV conversion)

    # --- 2. Calculate Molecular Properties ---
    # Convert masses to kg
    m_H_kg = m_H_amu * AMU_TO_KG
    m_Br_kg = m_Br_amu * AMU_TO_KG
    
    # Calculate reduced mass (mu)
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)
    
    # Calculate moment of inertia (I)
    I = mu * r0**2

    # --- 3. Calculate Spectroscopic Constants (in Joules) ---
    # Rotational constant B
    B_joules = (H_BAR**2) / (2 * I)
    
    # Centrifugal distortion constant D
    # Using the formula D = 4 * B^3 / (hbar^2 * (k/mu))
    D_joules = (4 * B_joules**3) / ((H_BAR**2) * (k / mu))
    
    # --- 4. Calculate Energy Shifts ---
    # Shift for J=0 -> J=1 is |ΔE| = 4 * D
    shift_0_to_1_joules = 4 * D_joules
    
    # Shift for J=1 -> J=2 is |ΔE| = 32 * D
    shift_1_to_2_joules = 32 * D_joules
    
    # --- 5. Convert to qeV ---
    # Conversion factor from Joules to qeV
    # shift_qeV = (shift_J / e) * 1e30
    JOULES_TO_QEV = (1 / E_CHARGE) * 1e30
    
    shift_0_to_1_qeV = shift_0_to_1_joules * JOULES_TO_QEV
    shift_1_to_2_qeV = shift_1_to_2_joules * JOULES_TO_QEV
    
    # --- 6. Print Results ---
    print("--- Intermediate Calculations ---")
    print(f"Reduced mass (μ): {mu:.4e} kg")
    print(f"Moment of inertia (I): {I:.4e} kg m^2")
    print(f"Rotational constant (B): {B_joules:.4e} J")
    print(f"Centrifugal distortion constant (D): {D_joules:.4e} J\n")

    print("--- Energy Shift for J=0 -> J=1 Transition ---")
    print(f"Equation for shift: |ΔE| = 4 * D")
    print(f"|ΔE| = 4 * {D_joules:.4e} J = {shift_0_to_1_joules:.4e} J")
    print(f"Energy Shift = {shift_0_to_1_qeV:.4e} qeV\n")

    print("--- Energy Shift for J=1 -> J=2 Transition ---")
    print(f"Equation for shift: |ΔE| = 32 * D")
    print(f"|ΔE| = 32 * {D_joules:.4e} J = {shift_1_to_2_joules:.4e} J")
    print(f"Energy Shift = {shift_1_to_2_qeV:.4e} qeV")

    return shift_0_to_1_qeV, shift_1_to_2_qeV

# Run the calculation and store the final answers
answer1, answer2 = calculate_centrifugal_distortion_shifts()
# The final answer format is not compatible with python code,
# but the calculation is done and you can see the result from the printed output.
# Example format: print(f'<<<{answer1:.4e}, {answer2:.4e}>>>')