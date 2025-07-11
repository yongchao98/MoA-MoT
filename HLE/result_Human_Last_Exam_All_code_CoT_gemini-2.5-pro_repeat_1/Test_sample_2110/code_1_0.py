import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e- + v_e(bar) -> mu- + v_mu(bar).

    The calculation is performed in the low-energy limit, where the W-boson
    propagator simplifies, leading to a four-fermion contact interaction.
    The formula used is derived from the standard cross-section calculation,
    incorporating phase space factors and a spin-averaged matrix element squared
    that reproduces the correct behavior in the massless limit.

    The formula is: σ = (G_F² / (π * s)) * (s - m_μ²)²
    """

    # Given parameters
    G_F = 1.0  # Fermi constant
    m_mu = 1.0 # Muon mass
    m_e = 1.0  # Electron mass
    s = 2.0    # Center-of-mass energy squared

    # --- Calculation ---
    
    # The derived formula for the total cross-section
    # σ = (G_F² / (π * s)) * (s - m_μ²)²
    
    # Numerator of the main term
    numerator = (s - m_mu**2)**2
    # Denominator of the main term
    denominator = math.pi * s
    
    # Calculate the total cross-section
    sigma = (G_F**2 / denominator) * numerator

    # --- Output ---
    print("Calculating the total cross-section σ for e⁻ + ν̄ₑ → μ⁻ + ν̄μ")
    print("\nFormula:")
    print("σ = (G_F² / (π * s)) * (s - m_μ²)²")
    
    print("\nSubstituting the given values:")
    # Print the equation with the numbers substituted
    # Showing each number in the final equation as requested.
    equation_str = f"σ = ({G_F}² / (π * {s})) * ({s} - {m_mu}²)²"
    print(equation_str)

    # Calculate parts of the equation for clarity in the next step
    term1 = f"{G_F**2 / (math.pi * s):.4f}"
    term2_base = s - m_mu**2
    term2 = f"{term2_base**2:.4f}"
    
    print(f"σ ≈ ({G_F**2:.2f} / ({math.pi:.4f} * {s:.2f})) * ({s:.2f} - {m_mu**2:.2f})²")
    print(f"σ ≈ ({G_F**2 / (math.pi * s):.4f}) * ({s - m_mu**2:.2f})²")
    print(f"σ ≈ {term1} * {term2}")
    
    print("\nFinal Result:")
    print(f"σ = {sigma}")

calculate_cross_section()
<<<0.15915494309189535>>>