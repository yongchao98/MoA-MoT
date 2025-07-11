import math

def calculate_proton_threshold_energy():
    """
    This function calculates the threshold energy for a proton reacting with a CMB
    photon to produce a Delta baryon, based on the provided physical model.
    """

    # 1. Define constants and parameters as per the problem statement.
    # Masses are in GeV.
    m_delta = 1.233  # Mass of Delta baryon
    # Proton mass, rounded to three decimal places.
    m_p = 0.938      # Mass of proton
    
    # Temperature in Kelvin.
    T = 2.73
    T_base = 2.7
    
    # Boltzmann's constant in GeV/K, rounded to three decimal places for the numeric part.
    k = 8.617e-14
    
    # Coefficients for the CMB photon energy calculation.
    C1_coeff = 2.71
    C2 = 1.0e-10  # Perturbative coefficient in GeV/K^2
    
    # 2. Calculate the mean CMB photon energy (E_gamma_mean).
    delta_T = T - T_base
    
    # Calculate the two components of the energy
    E_gamma_main_term = C1_coeff * k * T
    E_gamma_pert_term = C2 * (delta_T**2)
    
    E_gamma_mean = E_gamma_main_term + E_gamma_pert_term
    
    # 3. Calculate the proton threshold energy (E_p_thresh).
    # Using the formula for head-on collisions: E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)
    numerator = m_delta**2 - m_p**2
    denominator = 4 * E_gamma_mean
    
    E_p_thresh = numerator / denominator

    # 4. Output the results as requested.
    # Print the equation with the numerical values.
    print("Threshold Energy Calculation:")
    print(f"E_p_thresh = (m_delta^2 - m_p^2) / (4 * E_gamma_mean)")
    print(f"E_p_thresh = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma_mean:.4e})")
    
    # Print the final calculated value in the required format.
    print("\nFinal Answer:")
    final_answer_str = f"{E_p_thresh:.3e}"
    print(f"The average threshold energy for the proton is {final_answer_str} GeV.")

    return final_answer_str

# Execute the calculation and capture the final answer for the special format.
final_answer = calculate_proton_threshold_energy()

print(f"\n<<<{final_answer}>>>")
