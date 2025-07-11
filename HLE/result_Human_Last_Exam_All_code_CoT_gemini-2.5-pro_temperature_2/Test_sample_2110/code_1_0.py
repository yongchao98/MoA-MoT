import math

def calculate_cross_section():
    """
    Calculates the total cross-section for e- + ν̄e → μ- + ν̄μ
    in the low-energy limit.
    """
    # Given values
    G_F = 1.0
    m_mu = 1.0
    m_e = 1.0
    s = 2.0
    
    # Calculate squared masses for clarity
    m_mu_sq = m_mu**2
    m_e_sq = m_e**2
    
    print("Calculating the total cross-section σ using the formula:")
    print("σ = [ G_F² * (s - m_μ²)² * (2s² + s(m_e² + m_μ²) + 2m_e²m_μ²) ] / [ 12 * π * s³ ]\n")

    # Numerator calculation step-by-step
    term1_val = (s - m_mu_sq)**2
    term2_val = (2 * s**2 + s * (m_e_sq + m_mu_sq) + 2 * m_e_sq * m_mu_sq)
    
    numerator_val = G_F**2 * term1_val * term2_val

    # Denominator calculation
    denominator_val = 12 * math.pi * s**3
    
    # Final cross-section
    sigma = numerator_val / denominator_val
    
    # --- Printing the equation with values ---
    print("Substituting the given values:")
    print(f"G_F = {G_F}, m_e = {m_e}, m_μ = {m_mu}, s = {s}\n")

    # Building the string for the equation breakdown
    # Term 1: (s - m_μ²)²
    t1_str = f"({s} - {m_mu}²)²"
    # Term 2: (2s² + s(m_e² + m_μ²) + 2m_e²m_μ²)
    t2_str = f"(2*{s}² + {s}({m_e}² + {m_mu}²) + 2*{m_e}²*{m_mu}²)"
    
    # Top line of the fraction
    num_str_symbolic = f"{G_F}² * {t1_str} * {t2_str}"
    
    # Bottom line of the fraction
    den_str_symbolic = f"12 * π * {s}³"

    print(f"σ = [ {num_str_symbolic} ] / [ {den_str_symbolic} ]\n")
    
    # --- Now showing the evaluation of each part ---
    
    # Evaluate Term 1
    t1_eval = f"({s} - {m_mu_sq})² = {term1_val}"
    # Evaluate Term 2 parts
    t2p1_eval = 2*s**2
    t2p2_eval = s*(m_e_sq+m_mu_sq)
    t2p3_eval = 2*m_e_sq*m_mu_sq
    # Evaluate full Term 2
    t2_eval = f"(2*{s**2} + {s}({m_e_sq} + {m_mu_sq}) + 2*{m_e_sq}*{m_mu_sq}) = ({t2p1_eval} + {t2p2_eval} + {t2p3_eval}) = {term2_val}"

    # Evaluated top line of the fraction
    num_str_eval = f"{G_F**2} * {term1_val} * {term2_val}"
    # Evaluated bottom line of the fraction
    den_str_eval = f"12 * π * {s**3}"
    
    print("Evaluating the terms:")
    print(f"(s - m_μ²)² = {t1_eval}")
    print(f"(2s² + s(m_e² + m_μ²) + 2m_e²m_μ²) = {t2_eval}\n")

    print(f"σ = [ {num_str_eval} ] / [ {den_str_eval} ]")
    print(f"σ = {numerator_val} / {denominator_val:.4f}")
    print(f"σ = {numerator_val} / ({12*8}*π)")
    print(f"σ = {int(numerator_val)} / (96π) = {int(numerator_val/2)} / (48π) = 7 / (48π)\n")
    
    print("Final calculated cross-section:")
    print(f"σ ≈ {sigma}")


if __name__ == '__main__':
    calculate_cross_section()
    print("\n<<<0.04642055653158079>>>")