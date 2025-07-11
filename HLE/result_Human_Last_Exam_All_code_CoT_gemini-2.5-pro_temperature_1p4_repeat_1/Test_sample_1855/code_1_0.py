import math

def calculate_proton_threshold_energy():
    """
    Calculates the threshold energy for a proton reacting with a CMB photon to produce a Delta baryon.
    The calculation follows these steps:
    1. Define all physical constants with the specified precision.
    2. Calculate the average energy of a CMB photon at T=2.73 K, including the perturbative effect.
    3. Use the threshold energy formula derived from relativistic kinematics.
    4. Print the final equation with all numerical values substituted.
    5. Compute and print the final result in scientific notation, rounded to three decimal places.
    """

    # 1. Define constants and masses with specified precision.
    # Masses in GeV
    m_delta = 1.233
    m_proton = 0.938  # Rounded to three decimal places

    # Temperatures in Kelvin
    T = 2.73
    T_0 = 2.7

    # Physical constants
    # Per instructions, k's mantissa is rounded to 3 decimal places.
    k_boltzmann = 8.617e-14  # GeV/K
    
    # Coefficients for energy calculation
    c1_factor = 2.71
    c2_perturbative = 1.0e-10  # GeV/K^2

    # 2. Calculate the average energy of a CMB photon (E_gamma_bar)
    delta_T = T - T_0
    
    # E_gamma = (2.71 * k * T) + (C2 * delta_T^2)
    e_gamma_bar = (c1_factor * k_boltzmann * T) + (c2_perturbative * delta_T**2)

    # 3. Calculate the proton threshold energy (E_p_th)
    # E_p_th = (m_delta^2 - m_proton^2) / (4 * E_gamma_bar)
    numerator = m_delta**2 - m_proton**2
    denominator = 4 * e_gamma_bar
    e_p_th = numerator / denominator
    
    # 4. Print the final equation with numerical values
    print("The formula for the proton threshold energy is:")
    print("E_p,th = (m_delta^2 - m_proton^2) / (4 * E_gamma)")
    print("\nSubstituting the numerical values:")
    print(f"E_p,th = (({m_delta})^2 - ({m_proton})^2) / (4 * {e_gamma_bar})")

    # Evaluate the components of the equation
    print("\nWhere:")
    print(f"  m_delta^2 - m_proton^2 = {numerator:.6f} GeV^2")
    print(f"  4 * E_gamma = {denominator:.6e} GeV")
    
    # 5. Print the final result in the required format
    print("\nThe final threshold energy for the proton is:")
    final_answer_str = f"{e_p_th:.3e} GeV"
    print(final_answer_str)

    # Final answer in the specified format
    print(f"<<<{final_answer_str.split(' ')[0]}>>>")

# Run the calculation
calculate_proton_threshold_energy()