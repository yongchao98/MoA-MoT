import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e- + ν̄e -> μ- + ν̄μ.

    The formula used is:
    σ = (G_F^2 / (3*π*s)) * ((s - m_μ^2) / (s - m_e^2)) * (s^2 + (s/2)*(m_e^2 + m_μ^2) + m_e^2*m_μ^2)
    """
    # Given conditions
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0
    s = 2.0

    # Calculate intermediate terms for clarity in the final printout
    G_F_sq = G_F**2
    m_mu_sq = m_mu**2
    m_e_sq = m_e**2
    s_sq = s**2
    
    # Numerator of the ratio
    num_ratio = s - m_mu_sq
    # Denominator of the ratio
    den_ratio = s - m_e_sq
    # Ratio
    ratio = num_ratio / den_ratio

    # The term in the last bracket
    bracket_term = s_sq + (s / 2.0) * (m_e_sq + m_mu_sq) + m_e_sq * m_mu_sq

    # The prefactor
    prefactor = G_F_sq / (3 * math.pi * s)

    # Calculate the total cross-section
    sigma = prefactor * ratio * bracket_term

    # Print the step-by-step calculation with numbers
    print("The formula for the total cross-section σ is:")
    print("σ = (G_F^2 / (3*π*s)) * ((s - m_μ^2) / (s - m_e^2)) * (s^2 + s/2*(m_e^2 + m_μ^2) + m_e^2*m_μ^2)")
    print("\nSubstituting the given values:")
    print(f"G_F = {G_F}, m_μ = {m_mu}, m_e = {m_e}, s = {s}\n")
    
    # Output each number in the final equation as requested
    print("σ = ({g_f}^2 / (3*π*{s_val})) * (({s_val} - {m_mu}^2) / ({s_val} - {m_e}^2)) * ({s_val}^2 + {s_val}/2*({m_e}^2 + {m_mu}^2) + {m_e}^2*{m_mu}^2)".format(
        g_f=G_F, s_val=s, m_mu=m_mu, m_e=m_e
    ))
    
    # Show the evaluation of the parts
    print(f"σ = ({G_F_sq} / (3*π*{s})) * (({num_ratio}) / ({den_ratio})) * ({s_sq} + {s/2.0}*({m_e_sq} + {m_mu_sq}) + {m_e_sq*m_mu_sq})")
    print(f"σ = ({prefactor:.4f}) * ({ratio}) * ({bracket_term})")
    print(f"σ = {sigma}")

calculate_cross_section()
<<<0.3713678493133135>>>