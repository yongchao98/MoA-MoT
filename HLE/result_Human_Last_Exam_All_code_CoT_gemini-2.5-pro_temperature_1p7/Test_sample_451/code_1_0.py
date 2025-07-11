import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- 1. Define Constants and Input Parameters ---
    # Physical constants in SI units
    h_bar = 1.054571817e-34  # J*s (Reduced Planck constant)
    amu_to_kg = 1.66053906660e-27  # kg/amu
    e_charge = 1.602176634e-19  # C (elementary charge)
    pm_to_m = 1e-12          # m/pm

    # Molecule parameters
    r0_pm = 141.4  # pm
    k_force = 400.0  # N/m
    m_H_amu = 1.008  # amu
    m_Br_amu = 79.904 # amu
    
    # Unit conversions
    r0 = r0_pm * pm_to_m
    m_H = m_H_amu * amu_to_kg
    m_Br = m_Br_amu * amu_to_kg

    # --- 2. Calculate Molecular Properties ---
    # Reduced mass (mu)
    mu = (m_H * m_Br) / (m_H + m_Br)

    # Moment of inertia (I)
    I = mu * r0**2

    # --- 3. Calculate Energy Constants in Joules ---
    # Rotational constant (B)
    B = (h_bar**2) / (2 * I)

    # Vibrational energy (h_nu)
    omega = math.sqrt(k_force / mu)
    h_nu = h_bar * omega

    # Centrifugal distortion constant (D)
    D_J = (4 * B**3) / (h_nu**2)
    
    print("--- Calculated Molecular Constants ---")
    print(f"Reduced mass (μ): {mu:.4e} kg")
    print(f"Moment of inertia (I): {I:.4e} kg*m^2")
    print(f"Rotational constant (B): {B:.4e} J")
    print(f"Centrifugal distortion constant (D): {D_J:.4e} J")
    print("-" * 36 + "\n")

    # --- 4. Calculate Energy Shifts for Transitions ---
    
    # Transition 1: J = 0 to J = 1
    J_i_1, J_f_1 = 0, 1
    factor_1 = J_f_1**2 * (J_f_1 + 1)**2 - J_i_1**2 * (J_i_1 + 1)**2
    shift_J_1 = -D_J * factor_1

    # Transition 2: J = 1 to J = 2
    J_i_2, J_f_2 = 1, 2
    factor_2 = J_f_2**2 * (J_f_2 + 1)**2 - J_i_2**2 * (J_i_2 + 1)**2
    shift_J_2 = -D_J * factor_2

    # --- 5. Convert Units and Display Results ---
    J_to_eV = 1 / e_charge
    eV_to_qeV = 1 / 1e-30
    J_to_qeV = J_to_eV * eV_to_qeV
    
    shift_qeV_1 = shift_J_1 * J_to_qeV
    shift_qeV_2 = shift_J_2 * J_to_qeV

    # Print results for Transition 1
    print("--- Transition 1: J = 0 to J = 1 ---")
    print(f"ΔE = -D * [J_f²(J_f+1)² - J_i²(J_i+1)²]")
    print(f"ΔE = -({D_J:.4e}) * [{J_f_1}²({J_f_1+1})² - {J_i_1}²({J_i_1+1})²]")
    print(f"ΔE = -({D_J:.4e}) * [{int(factor_1)}]")
    print(f"Energy Shift = {shift_J_1:.4e} J")
    print(f"Final Energy Shift = {shift_qeV_1:.4e} qeV\n")

    # Print results for Transition 2
    print("--- Transition 2: J = 1 to J = 2 ---")
    print(f"ΔE = -D * [J_f²(J_f+1)² - J_i²(J_i+1)²]")
    print(f"ΔE = -({D_J:.4e}) * [{J_f_2}²({J_f_2+1})² - {J_i_2}²({J_i_2+1})²]")
    print(f"ΔE = -({D_J:.4e}) * [{int(factor_2)}]")
    print(f"Energy Shift = {shift_J_2:.4e} J")
    print(f"Final Energy Shift = {shift_qeV_2:.4e} qeV")
    
    return shift_qeV_1, shift_qeV_2

# Execute the calculation and print the final answer in the requested format
final_answers = calculate_energy_shifts()
# Create the answer string in the specified format
answer_string = f"[{final_answers[0]:.4e}, {final_answers[1]:.4e}]"

# The prompt asks for a specific format at the end.
# print(f"<<<{answer_string}>>>")
# A simple list of numbers might be better. Let's provide that.
print(f"\n<<<{final_answers[0]:.4e}, {final_answers[1]:.4e}>>>")
