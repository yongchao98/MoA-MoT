import math

def calculate_threshold_energy():
    """
    This function calculates the threshold energy for a proton reacting with a CMB photon
    to produce a Delta baryon, based on the provided physical model.
    """

    # Step 1: Define physical constants and given values.
    # Masses are given in GeV, and required constants are rounded to three decimal places.
    m_delta = 1.233  # Mass of Delta baryon in GeV
    m_p = 0.938      # Mass of proton in GeV (rounded from 938.272 MeV)
    k_B = 8.617e-14  # Boltzmann's constant in GeV/K (rounded from 8.617333e-5 eV/K)

    # Given parameters for the CMB photon energy calculation
    T = 2.73          # Current CMB temperature in K
    T0 = 2.7          # Base temperature in K
    C1_coeff = 2.71   # Coefficient for the linear energy term
    C2_coeff = 1.0e-10  # Coefficient for the perturbative energy term in GeV/K^2

    print("Defined constants and parameters:")
    print(f"  m_delta = {m_delta} GeV")
    print(f"  m_p = {m_p} GeV")
    print(f"  k_B = {k_B:.3e} GeV/K")
    print(f"  T = {T} K")
    print("-" * 35)

    # Step 2: Calculate the mean energy of the CMB photon (E_gamma).
    delta_T = T - T0
    E_gamma = (C1_coeff * k_B * T) + (C2_coeff * (delta_T**2))

    print("1. Calculate the mean CMB photon energy (E_gamma)")
    print("   Formula: E_gamma = (C1 * k * T) + (C2 * (T - T0)^2)")
    print(f"   E_gamma = ({C1_coeff} * {k_B:.3e} * {T}) + ({C2_coeff:.1e} * ({T} - {T0})^2)")
    
    # Calculate each part of the sum
    term1_E = C1_coeff * k_B * T
    term2_E = C2_coeff * delta_T**2
    
    print(f"   E_gamma = {term1_E:.6e} + {term2_E:.6e}")
    print(f"   E_gamma = {E_gamma:.6e} GeV\n")
    print("-" * 35)

    # Step 3: Calculate the threshold energy for the proton (E_p_thresh).
    # The exact relativistic formula for a head-on collision threshold is:
    # E_p_thresh = (m_final^2 - m_initial^2) / (4*E_photon) in the high-energy limit,
    # or more accurately, E_p_thresh = (m_delta^2 - m_p^2)/(4*E_gamma) + (m_p^2*E_gamma)/(m_delta^2-m_p^2)
    m_delta_sq = m_delta**2
    m_p_sq = m_p**2
    diff_sq_mass = m_delta_sq - m_p_sq
    
    term1_Ep = diff_sq_mass / (4 * E_gamma)
    term2_Ep = (m_p_sq * E_gamma) / diff_sq_mass
    E_p_thresh = term1_Ep + term2_Ep
    
    print("2. Calculate the proton's threshold energy (E_p)")
    print("   Formula: E_p = (m_delta^2 - m_p^2)/(4*E_gamma) + (m_p^2*E_gamma)/(m_delta^2-m_p^2)")
    print(f"   E_p = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma:.6e}) + ({m_p}^2 * {E_gamma:.6e}) / ({m_delta}^2 - {m_p}^2)")
    print(f"   E_p = ({m_delta_sq:.6f} - {m_p_sq:.6f}) / (4 * {E_gamma:.6e}) + ({m_p_sq:.6f} * {E_gamma:.6e}) / ({m_delta_sq:.6f} - {m_p_sq:.6f})")
    print(f"   E_p = {diff_sq_mass:.6f} / {4*E_gamma:.6e} + {m_p_sq * E_gamma:.6e} / {diff_sq_mass:.6f}")
    print(f"   E_p = {term1_Ep:.6e} + {term2_Ep:.6e}")
    print(f"   E_p = {E_p_thresh:.6e} GeV\n")
    print("-" * 35)
    
    # Step 4: Format and print the final answer.
    final_answer_str = f"{E_p_thresh:.3e}"
    
    print("Final Result:")
    print(f"The average threshold energy for the proton, rounded to three decimal places, is {final_answer_str} GeV.")

    return final_answer_str

# Execute the calculation and print the final answer in the required format
final_answer = calculate_threshold_energy()
print(f"<<<{final_answer}>>>")
