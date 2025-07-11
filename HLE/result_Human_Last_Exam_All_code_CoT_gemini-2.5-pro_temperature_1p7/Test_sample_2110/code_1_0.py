import math

def calculate_cross_section():
    """
    Calculates the total cross-section for e- + ν̄e → μ- + ν̄μ.

    The formula for the total cross-section (σ) is:
    σ = (G_F² / (π * s)) * (s - m_e²) * (s - m_μ²)

    This formula is derived from V-A theory in the low-energy limit.
    """
    # Given constants
    G_F = 1  # Fermi constant
    m_mu = 1 # Muon mass
    m_e = 1  # Electron mass
    s = 2    # Center-of-mass energy squared

    # Value of pi
    pi_val = math.pi

    # Step-by-step calculation
    # Term 1: G_F^2
    g_f_squared = G_F**2

    # Term 2: (s - m_e^2)
    s_minus_me_sq = s - m_e**2
    
    # Term 3: (s - m_mu^2)
    s_minus_mmu_sq = s - m_mu**2

    # Denominator: (pi * s)
    denominator = pi_val * s
    
    # Numerator of the main fraction
    numerator = g_f_squared * s_minus_me_sq * s_minus_mmu_sq

    # Calculate the total cross-section
    sigma = numerator / denominator

    # Output the equation with numerical values
    print("Calculating the total cross-section σ using the formula:")
    print("σ = (G_F² / (π * s)) * (s - m_e²) * (s - m_μ²)\n")
    print("Substituting the given values:")
    print(f"G_F = {G_F}")
    print(f"m_e = {m_e}")
    print(f"m_μ = {m_mu}")
    print(f"s = {s}")
    print(f"π ≈ {pi_val:.4f}\n")
    
    print("The equation becomes:")
    print(f"σ = ({g_f_squared} / ({pi_val:.4f} * {s})) * ({s} - {m_e**2}) * ({s} - {m_mu**2})")
    print(f"σ = ({g_f_squared} / {denominator:.4f}) * ({s_minus_me_sq}) * ({s_minus_mmu_sq})")
    print(f"σ = {numerator} / {denominator:.4f}\n")

    # Output the final result
    print("Final Result:")
    print(f"σ = {sigma}")

if __name__ == "__main__":
    calculate_cross_section()