import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e- + ν̄e → μ- + ν̄μ
    in the low-energy limit.
    """
    # Given parameters
    G_F = 1.0  # Fermi constant
    m_mu = 1.0 # Muon mass
    m_e = 1.0  # Electron mass (given, but not needed in the final formula)
    s = 2.0    # Center-of-mass energy squared

    # The theoretical formula for the total cross-section in the low-energy limit is:
    # σ = (G_F² * (s - m_μ²)²) / (π * s)

    # Calculate the numerator and denominator
    numerator = G_F**2 * (s - m_mu**2)**2
    denominator = math.pi * s

    # Calculate the total cross-section
    sigma = numerator / denominator

    # Print the final equation with all the numbers
    print("Cross-section calculation:")
    print("σ = (G_F² * (s - m_μ²)²) / (π * s)")
    print(f"σ = ({G_F:.1f}² * ({s:.1f} - {m_mu:.1f}²)²) / (π * {s:.1f})")
    print(f"σ = ({G_F**2:.1f} * ({s - m_mu**2:.1f})²) / ({math.pi:.5f} * {s:.1f})")
    print(f"σ = ({G_F**2:.1f} * ({(s - m_mu**2)**2:.1f})) / ({denominator:.5f})")
    print(f"σ = {numerator:.1f} / {denominator:.5f}")
    print(f"σ ≈ {sigma}")

if __name__ == '__main__':
    calculate_cross_section()
    # For the final answer extraction
    G_F = 1.0
    m_mu = 1.0
    s = 2.0
    sigma_val = (G_F**2 * (s - m_mu**2)**2) / (math.pi * s)
    # The final answer format requires only the numerical value
    # print(f'<<<{sigma_val}>>>')