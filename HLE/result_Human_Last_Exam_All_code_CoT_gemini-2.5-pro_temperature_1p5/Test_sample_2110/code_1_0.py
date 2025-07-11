import math

def calculate_cross_section(G_F, s, m_e, m_mu):
    """
    Calculates the total cross-section for the process e- + v_e(bar) -> mu- + v_mu(bar)
    in the low-energy limit.

    Args:
        G_F (float): Fermi constant.
        s (float): Center-of-mass energy squared.
        m_e (float): Electron mass.
        m_mu (float): Muon mass.
    
    Returns:
        float: The total cross-section sigma.
    """
    # Check for kinematic threshold
    if s <= m_mu**2 or s <= m_e**2:
        print("Energy s is below the threshold for the reaction.")
        return 0.0
        
    m_e_sq = m_e**2
    m_mu_sq = m_mu**2
    s_sq = s**2
    s_cub = s**3
    
    # The derived formula for the total cross-section is:
    # σ = (G_F² / (24*π*s³)) * (s - m_μ²)² * (2s² + s(m_e² + m_μ²) + 2*m_e²*m_μ²)
    
    term1_num = G_F**2
    term1_den = 24 * math.pi * s_cub
    
    term2 = (s - m_mu_sq)**2
    
    term3 = (2 * s_sq + s * (m_e_sq + m_mu_sq) + 2 * m_e_sq * m_mu_sq)
    
    sigma = (term1_num / term1_den) * term2 * term3
    
    # Print the equation with the given numbers
    print("Calculating the cross-section σ using the formula:")
    print("σ = (G_F² * (s - m_μ²)² / (24 * π * s³)) * (2s² + s(m_e² + m_μ²) + 2m_e²m_μ²)\n")
    print("Substituting the given values:")
    print(f"G_F = {G_F}")
    print(f"s = {s}")
    print(f"m_e = {m_e}")
    print(f"m_μ = {m_mu}\n")
    
    print("The final equation with numbers is:")
    print(f"σ = ({G_F}² * ({s} - {m_mu}²)² / (24 * π * {s}³)) * (2*{s}² + {s}({m_e}² + {m_mu}²) + 2*{m_e}²*{m_mu}²)")
    
    print(f"σ = ({term1_num} * {term2} / (24 * π * {s_cub})) * ({term3})")
    
    return sigma

# Given values
G_F = 1.0
s = 2.0
m_e = 1.0
m_mu = 1.0

# Calculate and print the result
cross_section = calculate_cross_section(G_F, s, m_e, m_mu)
print(f"\nFinal calculated cross-section σ = {cross_section}")
print(f"As a fraction of pi: 7 / (96 * π) ≈ {7/(96*math.pi)}")

# The problem asks to output just the final value.
# The calculation shows: 
# sigma = (1^2 * (2-1)^2 / (24*pi*2^3)) * (2*2^2 + 2*(1^2+1^2) + 2*1^2*1^2)
# sigma = (1 * 1 / (192*pi)) * (8 + 4 + 2)
# sigma = 14 / (192*pi) = 7 / (96*pi)

final_value = 7 / (96 * math.pi)
print(f"\n<<<__{final_value}__>>>")