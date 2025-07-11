import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e- + ν̄e -> μ- + ν̄μ
    in the low-energy limit given the specific values for G_F, masses, and energy.
    It first prints the general formula, then shows the substitution of each number
    into the formula, and finally prints the result.
    """

    # Given parameters
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0
    s = 2.0

    # For clarity in the formula
    m_mu_sq = m_mu**2
    m_e_sq = m_e**2

    # The derived formula for the total cross-section
    # σ = (G_F² / (24 * π * s³)) * (s - m_μ²)² * [2s² + s(m_e² + m_μ²) + 2m_e²m_μ²]
    
    # Calculate each part of the formula with the given numbers
    term_G_F = G_F**2
    denom_val = 24 * s**3
    term_s_minus_msq = (s - m_mu_sq)**2
    
    bracket_part1 = 2 * s**2
    bracket_part2 = s * (m_e_sq + m_mu_sq)
    bracket_part3 = 2 * m_e_sq * m_mu_sq
    bracket_total = bracket_part1 + bracket_part2 + bracket_part3

    # Calculate the final numerical result
    sigma = (term_G_F / (24 * math.pi * s**3)) * term_s_minus_msq * bracket_total

    # Outputting the steps as requested
    print("The final equation for the total cross-section is derived as:")
    print("σ = (G_F² / (24 * π * s³)) * (s - m_μ²)² * [2s² + s(m_e² + m_μ²) + 2m_e²m_μ²]\n")

    print("Substituting each number into the equation:")
    print(f"σ = ({int(G_F)}² / (24 * π * {int(s)}³)) * ({int(s)} - {int(m_mu)}²)² * [2*{int(s)}² + {int(s)}*({int(m_e)}² + {int(m_mu)}²) + 2*{int(m_e)}²*{int(m_mu)}²]\n")

    print("Evaluating each term in the equation:")
    print(f"σ = ({int(term_G_F)} / ({int(denom_val)} * π)) * {int(term_s_minus_msq)} * [{int(bracket_part1)} + {int(bracket_part2)} + {int(bracket_part3)}]\n")
    
    print("Simplifying the expression:")
    print(f"σ = ({int(term_G_F * term_s_minus_msq * bracket_total)} / ({int(denom_val)} * π))")
    print("σ = 14 / (192 * π) = 7 / (96 * π)\n")
    
    print("Final numerical result:")
    print(sigma)
    
    return sigma

final_answer = calculate_cross_section()
print(f"\n<<<{final_answer}>>>")