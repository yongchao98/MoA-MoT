import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e- + v_e(bar) -> mu- + v_mu(bar)
    in the low-energy limit.
    """
    # Given parameters
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0
    s = 2.0

    # Squares of masses
    m_mu_sq = m_mu**2
    m_e_sq = m_e**2

    print("The formula for the total cross-section (σ) is:")
    print("σ = (G_F² * (s - m_μ²)²) / (12 * π * s³) * [2*s² + s*(m_e² + m_μ²) + 2*m_e²*m_μ²]")
    print("\nPlugging in the values:")
    print(f"G_F = {G_F}, m_e = {m_e}, m_μ = {m_mu}, s = {s}")
    print(f"m_e² = {m_e_sq}, m_μ² = {m_mu_sq}")
    
    # Step-by-step evaluation of the formula with the given numbers
    print("\nEquation with numbers:")
    # The user wants to see each number in the equation.
    print(f"σ = ({G_F**2:.1f} * ({s:.1f} - {m_mu_sq:.1f})²) / (12 * π * {s:.1f}³) * [2*({s:.1f}²) + {s:.1f}*({m_e_sq:.1f} + {m_mu_sq:.1f}) + 2*{m_e_sq:.1f}*{m_mu_sq:.1f}]")

    # Calculate intermediate terms to show the simplification
    s_minus_m_mu_sq_sq = (s - m_mu_sq)**2
    s_cubed = s**3
    s_squared = s**2
    m_sum_sq = m_e_sq + m_mu_sq
    m_prod_sq = m_e_sq * m_mu_sq
    
    print("\nSimplifying the terms:")
    print(f"σ = ({G_F**2 * s_minus_m_mu_sq_sq:.1f}) / (12 * π * {s_cubed:.1f}) * [2*{s_squared:.1f} + {s:.1f}*{m_sum_sq:.1f} + 2*{m_prod_sq:.1f}]")

    term_in_bracket_1 = 2 * s_squared
    term_in_bracket_2 = s * m_sum_sq
    term_in_bracket_3 = 2 * m_prod_sq
    bracket_total = term_in_bracket_1 + term_in_bracket_2 + term_in_bracket_3
    
    print(f"σ = {G_F**2 * s_minus_m_mu_sq_sq:.1f} / ({12 * s_cubed:.1f} * π) * [{term_in_bracket_1:.1f} + {term_in_bracket_2:.1f} + {term_in_bracket_3:.1f}]")
    print(f"σ = {G_F**2 * s_minus_m_mu_sq_sq:.1f} / ({12 * s_cubed:.1f} * π) * {bracket_total:.1f}")
    
    numerator = G_F**2 * s_minus_m_mu_sq_sq * bracket_total
    denominator_coeff = 12 * s_cubed
    
    print(f"σ = {numerator:.1f} / ({denominator_coeff:.1f} * π)")
    
    # Simplify the fraction 14/96 to 7/48
    import fractions
    f = fractions.Fraction(int(numerator), int(denominator_coeff))
    print(f"σ = {f.numerator} / ({f.denominator} * π)")

    # Final numerical calculation
    sigma = (G_F**2 * (s - m_mu_sq)**2 / (12 * math.pi * s**3)) * (2*s**2 + s*(m_e_sq + m_mu_sq) + 2*m_e_sq*m_mu_sq)

    print("\nFinal calculated value:")
    print(f"σ ≈ {sigma}")

calculate_cross_section()
<<<0.04642118615872341>>>