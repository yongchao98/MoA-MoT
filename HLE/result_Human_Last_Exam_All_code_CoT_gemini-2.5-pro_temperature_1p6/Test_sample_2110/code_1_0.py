import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e⁻ + ν̄ₑ → μ⁻ + ν̄μ
    in the low-energy limit.
    """
    # Given values
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0
    s = 2.0
    pi = math.pi

    # The derived formula for the total cross-section is:
    # σ = (G_F² / (48πs³)) * (s - m_μ²)² * [4s² + 2s(m_e² + m_μ²) + 4m_e²m_μ²]

    # Calculate each part of the formula with the given values
    # for the purpose of demonstrating the final equation with numbers
    gf_sq = G_F**2
    s_cubed = s**3
    coeff_val = 48 * pi * s_cubed
    
    muon_term_val = (s - m_mu**2)**2

    s_sq_term_val = 4 * s**2
    mass_sum_val = m_e**2 + m_mu**2
    mixed_term_val = 2 * s * mass_sum_val
    mass_prod_val = 4 * m_e**2 * m_mu**2
    bracket_term_val = s_sq_term_val + mixed_term_val + mass_prod_val

    # Print the equation with substituted values
    print("The total cross-section σ is calculated using the formula:")
    print("σ = (G_F² / (48 * π * s³)) * (s - m_μ²)² * (4s² + 2s(m_e² + m_μ²) + 4m_e²m_μ²)\n")
    print("Substituting the given values:")
    print(f"G_F = {G_F}, m_e = {m_e}, m_μ = {m_mu}, s = {s}\n")
    
    print("The equation becomes:")
    # Showing each number in the final equation
    print(f"σ = ({G_F:.1f}² / (48 * π * {s:.1f}³)) * ({s:.1f} - {m_mu:.1f}²)² * (4*{s:.1f}² + 2*{s:.1f}({m_e:.1f}² + {m_mu:.1f}²) + 4*{m_e:.1f}²*{m_mu:.1f}²)")
    
    # Show the calculation of the terms
    print(f"σ = ({gf_sq:.1f} / ({coeff_val:.2f})) * ({muon_term_val:.1f}) * ({s_sq_term_val:.1f} + {mixed_term_val:.1f} + {mass_prod_val:.1f})")
    print(f"σ = ({gf_sq / coeff_val:.6f}) * ({muon_term_val:.1f}) * ({bracket_term_val:.1f})")

    # Calculate the final result
    sigma = (gf_sq / coeff_val) * muon_term_val * bracket_term_val
    
    print(f"\nThe calculated total cross-section is:")
    print(f"σ = {sigma}")

calculate_cross_section()
<<<0.023183181813977533>>>