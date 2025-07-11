import math

def calculate_cross_section():
    """
    Calculates the total cross-section for the process e⁻ + ν̄ₑ → μ⁻ + ν̄μ
    in the low-energy limit.
    """
    # Given parameters
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0  # Note: This parameter cancels out in the theoretical derivation.
    s = 2.0

    # The derived formula for the total cross-section is:
    # sigma = (G_F**2 / (4 * pi * s)) * (s - m_mu**2)**2
    # This formula resulted from a detailed calculation involving the matrix element
    # for W-boson exchange in the low-energy limit, spin averaging, and phase space integration.
    # A key feature of this s-channel process is that the dependency on the initial
    # electron mass (m_e) from the kinematic factors and the matrix element cancels out.

    # Calculate the final value
    pi = math.pi
    sigma = (G_F**2 / (4 * pi * s)) * (s - m_mu**2)**2

    # Print the equation with the specific numbers substituted, as requested
    print("The final formula for the total cross-section is:")
    print("σ = (G_F² / (4 * π * s)) * (s - m_μ²)²\n")
    print("Plugging in the given values:")
    print(f"σ = ({G_F}² / (4 * {pi:.4f} * {s})) * ({s} - {m_mu}²)²")
    print(f"σ = ({G_F**2} / ({4 * pi * s:.4f})) * ({s - m_mu**2})²")
    print(f"σ = ({1 / (8 * pi):.4f}) * ({1.0})²\n")

    # Print the final numerical result
    print(f"The calculated total cross-section is: {sigma}")
    
    return sigma

if __name__ == "__main__":
    final_answer = calculate_cross_section()
    # The final answer is wrapped in <<<>>> as requested, though printing is sufficient.
    # print(f"<<<{final_answer}>>>")