import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the electroweak interaction
    e⁻ + ν̄ₑ → μ⁻ + ν̄μ in the low-energy limit.

    The specified process is non-standard. We use the formula for the analogous,
    physically allowed process ν_μ + e⁻ → μ⁻ + ν_e, which is given by:
    σ = G_F² * (s - m_μ²)² / (π * s)

    This formula is valid under the low-energy approximation where the interaction
    is treated as a contact four-fermion interaction.

    Args are provided by the user in the prompt.
    G_F: Fermi constant
    m_mu: Muon mass
    m_e: Electron mass (note: this formula assumes m_e << sqrt(s))
    s: Center-of-mass energy squared
    """

    # Given values
    G_F = 1
    m_mu = 1
    m_e = 1 # Not used in the final formula, but given in the problem
    s = 2
    pi = math.pi

    # The formula for the total cross-section
    # σ = G_F² * (s - m_μ²)² / (π * s)
    sigma = (G_F**2 * (s - m_mu**2)**2) / (pi * s)

    # Output the calculation steps as requested by the user
    print(f"Calculating the total cross-section σ using the formula:")
    print(f"σ = G_F^2 * (s - m_μ^2)^2 / (π * s)\n")
    print(f"Given values:")
    print(f"G_F = {G_F}")
    print(f"m_μ = {m_mu}")
    print(f"s = {s}\n")
    print(f"Plugging the values into the equation:")
    # The prompt requests that each number in the final equation be output.
    print(f"σ = {G_F}^2 * ({s} - {m_mu}^2)^2 / (π * {s})")
    
    # Calculate and display the final result
    numerator = G_F**2 * (s - m_mu**2)**2
    denominator = pi * s
    print(f"σ = {numerator} / {denominator}")
    print(f"σ = {sigma}")

calculate_cross_section()