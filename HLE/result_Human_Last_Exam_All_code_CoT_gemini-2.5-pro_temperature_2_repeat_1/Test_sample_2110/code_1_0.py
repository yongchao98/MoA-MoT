import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e- + nubar_e -> mu- + nubar_mu
    in the low-energy limit.
    """
    # Given parameters
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0
    s = 2.0
    pi = math.pi

    # The formula for the total cross-section (sigma)
    # sigma = (G_F**2 / (pi * s)) * (s - m_e**2) * (s - m_mu**2)

    # Perform the calculation
    numerator = G_F**2 * (s - m_e**2) * (s - m_mu**2)
    denominator = pi * s
    sigma = numerator / denominator

    # Print the equation with variables
    print("The formula for the total cross-section is:")
    print("σ = (G_F² / (π * s)) * (s - m_e²) * (s - m_μ²)\n")

    # Print the equation with the plugged-in values
    print("Plugging in the given values:")
    print(f"σ = ({G_F}² / (π * {s})) * ({s} - {m_e}²) * ({s} - {m_mu}²)")
    
    # Calculate the intermediate steps for clarity
    term1 = G_F**2
    term2 = pi * s
    term3 = s - m_e**2
    term4 = s - m_mu**2
    
    print(f"σ = ({term1} / {term2:.4f}) * ({term3}) * ({term4})")
    print(f"σ = ({term1 / term2:.4f}) * ({term3 * term4})")
    
    # Print the final result
    print("\nThe final calculated total cross-section is:")
    print(f"σ = {sigma}")
    print(f"<<<{sigma}>>>")

calculate_cross_section()