import math

def calculate_cross_section():
    """
    Calculates the total cross-section for e- + v_e_bar -> mu- + v_mu_bar.
    """
    # Given parameters
    G_F = 1.0
    s = 2.0
    m_e = 1.0
    m_mu = 1.0

    # The formula for the total cross-section is:
    # sigma = (G_F^2 / (96 * pi * s^3)) * (s - m_mu^2)^2 * [2*s^2 + s*(m_e^2 + m_mu^2) + 2*m_e^2*m_mu^2]
    
    # Check if the process is physically possible (s > m_mu^2)
    if s <= m_mu**2:
        print("The center-of-mass energy squared (s) is not sufficient to produce a muon.")
        print("Process is not possible under these conditions (s <= m_mu^2).")
        return

    # Calculate each part of the formula with the given numbers
    term1 = G_F**2
    term2_numerator = (s - m_mu**2)**2
    term2_denominator = (96 * math.pi * s**3)
    
    bracket_term_1 = 2 * s**2
    bracket_term_2 = s * (m_e**2 + m_mu**2)
    bracket_term_3 = 2 * m_e**2 * m_mu**2
    bracket_term_total = bracket_term_1 + bracket_term_2 + bracket_term_3
    
    # Print the equation with all the numbers substituted in
    print("Calculating the total cross-section sigma (σ) using the formula:")
    print("σ = (G_F² / (96 * π * s³)) * (s - m_μ²)² * [2s² + s(mₑ² + m_μ²) + 2mₑ²m_μ²]\n")

    print(f"Substituting the values: G_F = {G_F}, s = {s}, mₑ = {m_e}, m_μ = {m_mu}\n")
    
    print("σ = ({}² / (96 * π * {}³)) * ({} - {}²)² * [2*{}² + {}*({}² + {}²) + 2*{}²*{}²]".format(
        G_F, s, s, m_mu, s, s, m_e, m_mu, m_e, m_mu))

    # Calculate the numerical value of the parts
    calc_part1 = G_F**2 / (96 * math.pi * s**3)
    calc_part2 = (s - m_mu**2)**2
    calc_part3 = (2 * s**2 + s * (m_e**2 + m_mu**2) + 2 * m_e**2 * m_mu**2)
    
    print(f"σ = ({calc_part1:.4f}) * ({calc_part2}) * ({calc_part3})")
    
    # Final result calculation
    sigma = calc_part1 * calc_part2 * calc_part3
    
    print("\nSimplifying the expression:")
    print("Numerator = {}² * ({} - {})² * [{} + {} + {}] = {} * {} * {} = {}".format(
        G_F, s, m_mu**2, bracket_term_1, bracket_term_2, bracket_term_3, term1, term2_numerator, bracket_term_total, term1 * term2_numerator * bracket_term_total))
    print("Denominator = 96 * π * {}³ = 96 * π * {} = {}".format(
        s, s**3, 96*s**3))
    print(f"σ = {term1 * term2_numerator * bracket_term_total} / ({96*s**3} * π) = {int(term1 * term2_numerator * bracket_term_total)} / {int(96*s**3)}π")
    print(f"σ = {7} / {384}π")


    print("\nFinal Result:")
    print(f"σ = {sigma}")


calculate_cross_section()
<<<0.005786591325145781>>>