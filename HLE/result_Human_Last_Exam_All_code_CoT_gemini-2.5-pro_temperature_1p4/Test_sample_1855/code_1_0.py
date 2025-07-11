# -*- coding: utf-8 -*-

def calculate_proton_threshold_energy():
    """
    This script calculates the threshold energy for a proton colliding with a CMB photon
    to produce a Delta baryon, including a perturbative effect on the CMB energy.
    """
    
    # Step 1: Define physical constants and given parameters.
    # Masses are in GeV, Boltzmann's constant is in GeV/K.
    # Physical constants are rounded to three decimal places as requested.
    m_delta = 1.233  # Mass of Delta baryon in GeV
    m_p = 0.938      # Mass of proton in GeV (rounded to 3 decimal places)
    k_B = 8.617e-14  # Boltzmann's constant in GeV/K (rounded to 3 decimal places)

    # Parameters for the CMB photon energy calculation
    T = 2.73         # CMB temperature in Kelvin
    T0 = 2.7         # Base temperature for perturbation in Kelvin
    c1_coeff = 2.71  # Proportionality coefficient for the linear energy term
    c2_coeff = 1.0e-10 # Coefficient for the perturbative energy term in GeV/K^2
    
    # Step 2: Calculate the mean energy of a CMB photon (E_gamma).
    delta_T = T - T0
    E_gamma = (c1_coeff * k_B * T) + (c2_coeff * delta_T**2)
    
    # Step 3: Calculate the proton's threshold energy (E_p_threshold).
    # The formula is derived from relativistic kinematics for a head-on collision:
    # E_p_th = (m_delta^2 - m_p^2) / (4 * E_gamma)
    E_p_threshold = (m_delta**2 - m_p**2) / (4 * E_gamma)
    
    # Step 4: Output the results as requested.
    # First, print the equation with the numerical values plugged in.
    print("The threshold energy for the proton is calculated with the following equation:")
    # This line explicitly shows the numbers used in the final equation.
    print(f"E_p_th = (({m_delta})^2 - ({m_p})^2) / (4 * {E_gamma:.6e})\n")
    
    # Format the final answer in scientific notation with three decimal places.
    final_answer_formatted = f"{E_p_threshold:.3e}"
    
    print("The final result for the proton threshold energy is:")
    print(f"{final_answer_formatted} GeV")
    
    # The final answer in the specified format for parsing.
    print(f"<<<{final_answer_formatted}>>>")

# Execute the calculation
calculate_proton_threshold_energy()