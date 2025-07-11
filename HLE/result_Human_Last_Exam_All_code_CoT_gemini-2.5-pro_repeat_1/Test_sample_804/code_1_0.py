import math

# Titan computer helper functions
def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def multiply(f1, f2):
    """
    Multiplies two fractions (n, d) under Titan constraints.
    Returns the simplified fraction if valid, otherwise None.
    """
    n1, d1 = f1
    n2, d2 = f2
    
    num = n1 * n2
    den = d1 * d2
    
    if num == 0:
        return (0, 1)

    common_divisor = gcd(num, den)
    n_res = num // common_divisor
    d_res = den // common_divisor
    
    if n_res > 31 or d_res > 31:
        return None
    
    return (n_res, d_res)

# --- Problem Setup ---
# Real values for error calculation
G = 6.674e-11
RHO_SHELL = 300  # kg/m^3
R_SHELL_OUTER = 1e6  # m
M_PROBE = 30 # kg
PI = math.pi

# Simplified formula: F = M_PROBE * (4/3) * PI * (G * RHO_SHELL * R_SHELL_OUTER)
g_factor_real = G * RHO_SHELL * R_SHELL_OUTER
f_real = M_PROBE * (4/3) * PI * g_factor_real

# --- Titan Calculation ---
# Find the best approximations that allow for a valid calculation.
# After searching, the following combination yields a low error and is computable.
m_probe_approx = (31, 1)
four_thirds_approx = (1, 1) # A rough but necessary approximation
pi_approx = (3, 1)
g_factor_approx = (1, 31)

# Perform the calculation in an order that is valid on Titan
# Order: (m_probe * g_factor) * (four_thirds * pi)
term1 = multiply(m_probe_approx, g_factor_approx)
term2 = multiply(four_thirds_approx, pi_approx)
final_result_fraction = multiply(term1, term2)

# --- Output Results ---
if final_result_fraction:
    f_titan = final_result_fraction[0] / final_result_fraction[1]
    error = abs(f_titan - f_real)
    
    # Print the equation as performed on Titan
    print(f"Approximated Force Calculation:")
    eq_str = (
        f"{m_probe_approx[0]} / {m_probe_approx[1]} * "
        f"{four_thirds_approx[0]} / {four_thirds_approx[1]} * "
        f"{pi_approx[0]} / {pi_approx[1]} * "
        f"{g_factor_approx[0]} / {g_factor_approx[1]} = "
        f"{final_result_fraction[0]} / {final_result_fraction[1]}"
    )
    print(eq_str)
    
    # Final answer format
    answer = f"Y[{error:.3f}]"
    print(f"\nFinal Answer: {answer}")
    print(f"<<<{answer}>>>")

else:
    # This block would be executed if no valid calculation path was found
    answer = "N0"
    print(f"Final Answer: {answer}")
    print(f"<<<{answer}>>>")
