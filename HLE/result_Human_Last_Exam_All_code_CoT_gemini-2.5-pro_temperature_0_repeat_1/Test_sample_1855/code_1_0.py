def calculate_threshold_energy():
    """
    Calculates the threshold energy for a proton reacting with a CMB photon
    to produce a Delta baryon.
    """
    # Step 1: Define all physical constants and given values.
    # Masses are in GeV, temperatures in Kelvin, and k in GeV/K.
    # Constants are rounded to three decimal places as specified.
    m_delta = 1.233      # Mass of Delta baryon in GeV
    m_p = 0.938          # Mass of proton in GeV
    k = 8.617e-14        # Boltzmann's constant in GeV/K
    
    # Given experimental parameters
    T = 2.73             # CMB temperature in K
    T_ref = 2.7          # Reference temperature for perturbation in K
    prop_coeff = 2.71    # Proportionality coefficient for the main energy term
    pert_coeff = 1.0e-10 # Perturbative coefficient in GeV/K^2

    # Step 2: Calculate the mean energy of a CMB photon (E_gamma).
    # The formula includes a primary term and a perturbative correction.
    # E_gamma = (2.71 * k * T) + C * (T - T_ref)^2
    delta_T = T - T_ref
    main_energy_term = prop_coeff * k * T
    perturbative_term = pert_coeff * (delta_T**2)
    E_gamma = main_energy_term + perturbative_term

    print("1. Calculating the mean CMB photon energy (E_gamma):")
    print(f"E_gamma = ({prop_coeff} * {k:.3e} GeV/K * {T} K) + ({pert_coeff:.1e} GeV/K^2 * ({T} K - {T_ref} K)^2)")
    print(f"E_gamma = {main_energy_term} GeV + {perturbative_term} GeV")
    print(f"E_gamma = {E_gamma} GeV\n")

    # Step 3: Calculate the proton's threshold energy (E_p_th).
    # The formula is derived from relativistic kinematics: E_p_th = (m_delta^2 - m_p^2) / (4 * E_gamma)
    numerator = m_delta**2 - m_p**2
    denominator = 4 * E_gamma
    E_p_th = numerator / denominator

    print("2. Calculating the proton threshold energy (E_p_th):")
    print(f"E_p_th = (m_delta^2 - m_p^2) / (4 * E_gamma)")
    print(f"E_p_th = ({m_delta}^2 GeV^2 - {m_p}^2 GeV^2) / (4 * {E_gamma} GeV)")
    print(f"E_p_th = ({m_delta**2} - {m_p**2}) GeV^2 / ({denominator} GeV)")
    print(f"E_p_th = {numerator} GeV^2 / {denominator} GeV")
    print(f"E_p_th = {E_p_th} GeV\n")

    # Step 4: Print the final answer in the required format.
    print("Final Answer:")
    print(f"The average threshold energy for the proton is {E_p_th:.3e} GeV.")

if __name__ == "__main__":
    calculate_threshold_energy()