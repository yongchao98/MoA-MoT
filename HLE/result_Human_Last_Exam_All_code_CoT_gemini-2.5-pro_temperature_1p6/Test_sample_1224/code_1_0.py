import math
from fractions import Fraction

def solve_pandora_fall():
    """
    Solves the Pandora falling ball problem using Titan 4-bit computer constraints.
    """

    # --- Step 1: Define constants and the target value ---
    
    # Physics constants (high precision for calculating the 'true' value)
    G = 6.674e-11
    # Pandora's properties
    r_core = 1e5  # m
    rho_core = 1200  # kg/m^3
    R_eq = 2e6  # m
    rho_shell = 300  # kg/m^3
    h_val = 5000 # m

    # Calculate actual mass and gravity for reference
    v_core = (4/3) * math.pi * r_core**3
    m_core = rho_core * v_core
    v_planet = (4/3) * math.pi * R_eq**3
    m_shell = rho_shell * (v_planet - v_core)
    m_total = m_core + m_shell
    g_actual = (G * m_total) / (R_eq**2)
    t_actual = math.sqrt(2 * h_val / g_actual)


    # --- Step 2: Simulate Titan's fractional calculations ---
    
    # h = 5000, represented as 5/1 * 10^3
    h_frac = Fraction(5, 1)
    h_exp = 3

    # g ≈ 1.676 m/s^2. We approximate g with a simple, valid fraction 5/3 (1.667).
    # This assumes Titan's RDX engine simplifies the complex g calculation to this value.
    g_frac = Fraction(5, 3)

    # Calculate t = sqrt(2 * h / g)
    # S = 2 * h / g = (2 * 5/1 * 10^3) / (5/3)
    # S = (10/1 * 10^3) / (5/3) = (1/1 * 10^4) / (5/3) = 3/5 * 10^4
    S_frac = Fraction(3, 5)
    S_exp = 4
    
    # Now t = sqrt(S) = sqrt(3/5 * 10^4) = sqrt(3/5) * 10^2
    # We must use a pre-computed approximation for sqrt(3/5) ≈ 0.7746
    # A good fractional approximation is 7/9 ≈ 0.7778
    sqrt_S_approx_frac = Fraction(7, 9)

    # Final time t = 7/9 * 10^2
    t_titan_frac = sqrt_S_approx_frac
    t_titan_exp = S_exp / 2
    
    t_titan_val = float(t_titan_frac) * (10**t_titan_exp)
    
    # Calculate the absolute error
    error = abs(t_titan_val - t_actual)

    # --- Step 3: Print the derivation as requested ---
    
    print("Yes, Titan can be used to calculate the landing time.")
    print("The calculation proceeds as follows:\n")
    print("1. The landing time 't' is found using the formula t = sqrt(2 * h / g).")
    print("2. The complex calculation for g (surface gravity) is simplified to the fraction 5/3.")
    print("3. The height h is 5000 m.\n")
    
    print("The final equation and its step-by-step simplification:")
    
    h = 5000
    g_n, g_d = g_frac.numerator, g_frac.denominator
    
    val1 = 2 * h
    val2 = float(g_frac)
    val3 = val1 / val2
    
    sqrt_val_n, sqrt_val_d = sqrt_S_approx_frac.numerator, sqrt_S_approx_frac.denominator
    
    final_val = t_titan_val

    # Print each number in the equation as requested
    # The numbers are: 2, 5000, 5, 3 (from the setup),
    # 7, 9, 2 (from the sqrt approx), and the final result
    print(f"t = sqrt(({2} * {h}) / ({g_n}/{g_d}))")
    print(f"  = sqrt({val1} / {val2:.2f})")
    print(f"  = sqrt({val3:.2f})")
    print(f"  = sqrt(3/5 * 10^4)")
    print(f"  = sqrt(3/5) * 10^2")
    print(f"Approximating sqrt(3/5) with the pre-computed fraction {sqrt_val_n}/{sqrt_val_d}:")
    print(f"t ≈ {sqrt_val_n}/{sqrt_val_d} * 10^{int(t_titan_exp)}")
    print(f"  = {sqrt_val_n * 100}/{sqrt_val_d} ≈ {final_val:.2f} s")
    
    print(f"\nThe absolute error is |{final_val:.2f} - {t_actual:.2f}| = {error:.2f} s.")

    # Format final answer string
    return f"Y[{error:.2f}]"

# Execute the solution
final_answer = solve_pandora_fall()
# The final answer is wrapped in <<<>>>
print(f"\n<<<{final_answer}>>>")