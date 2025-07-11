import math

# High-precision constants for true value calculation
G_TRUE = 6.67430e-11
PI_TRUE = math.pi
R_CORE_TRUE = 1.0e5  # m
R_TOTAL_TRUE = 2.0e6  # m
RHO_CORE_TRUE = 1200  # kg/m^3
RHO_SHELL_TRUE = 300  # kg/m^3

def calculate_true_v_escape():
    """Calculates the true escape velocity using high-precision float values."""
    v_core = (4/3) * PI_TRUE * (R_CORE_TRUE**3)
    m_core = RHO_CORE_TRUE * v_core

    v_total = (4/3) * PI_TRUE * (R_TOTAL_TRUE**3)
    v_shell = v_total - v_core
    m_shell = RHO_SHELL_TRUE * v_shell
    
    m_total = m_core + m_shell
    
    v_escape_sq = (2 * G_TRUE * m_total) / R_TOTAL_TRUE
    return math.sqrt(v_escape_sq)

def calculate_titan_v_escape():
    """
    Simulates the escape velocity calculation using Titan's 4-bit fractional arithmetic.
    This involves simplification and approximation.
    """
    # Simplified formula: v_e^2 approx 32 * pi * G * 10^14
    
    # Step 1: Choose best fractional approximations for pi and G
    # These values are chosen to allow cancellation and minimize error.
    pi_num, pi_den = 13, 4  # pi approx 3.25
    g_num, g_den = 13, 2  # G approx 6.5e-11

    # Step 2: Calculate the intermediate term 32 * pi * G
    # On Titan, this would involve ordered multiplication and cancellation, e.g., (32/4)*13*(13/2)
    # The result of this sub-calculation is 676.
    term_32_pi_g = 32 * (pi_num / pi_den) * (g_num / g_den)

    # Step 3: Calculate v_escape^2
    # v_e^2 approx 676 * 10^-11 * 10^14 = 676 * 10^3
    v_escape_sq_titan = term_32_pi_g * 1e3

    # Step 4: Approximate the square root with a 4-bit fraction
    # We need to find the best a/b for sqrt(67.6) ~ 8.222
    # Let's search all possible 4-bit fractions (1-15 for num/den)
    
    target = math.sqrt(v_escape_sq_titan) / 100 # Target is sqrt(67.6)
    best_frac_val = 0
    min_err = float('inf')
    best_num, best_den = 0, 0

    for a in range(1, 16):
        for b in range(1, 16):
            val = a / b
            if abs(val - target) < min_err:
                min_err = abs(val - target)
                best_frac_val = val
                best_num, best_den = a, b
    
    # The final approximated escape velocity
    v_escape_titan = best_frac_val * 100
    
    # Final display
    print("Titan Calculation for Pandora's Escape Velocity")
    print("-" * 50)
    print("Simplified equation: v_e^2 = 32 * pi * G * 10^14")
    print(f"Approximations: pi = {pi_num}/{pi_den}, G = {g_num}/{g_den}e-11")
    
    # To satisfy the "output each number" requirement, let's build the equation string
    # from the derived values.
    # v_e = sqrt(32 * pi * G * 1e14)
    # The final step is approximating sqrt with the best fraction
    print(f"\nFinal calculation step:")
    print(f"v_e = Sqrt_Approx( Sqrt( 32 * ({pi_num}/{pi_den}) * ({g_num}/{g_den}e-11) * 1e14 ) )")
    print(f"v_e = Sqrt_Approx( Sqrt( {v_escape_sq_titan:.2f} ) )")
    print(f"v_e = Sqrt_Approx( {math.sqrt(v_escape_sq_titan):.2f} )")
    print(f"Best single fraction approximation for sqrt(67.6) is {best_num}/{best_den}")
    print(f"Final approximated result: v_e = ({best_num}/{best_den}) * 100 = {v_escape_titan:.2f} m/s")

    return v_escape_titan

v_true = calculate_true_v_escape()
v_titan = calculate_titan_v_escape()

abs_error = abs(v_true - v_titan)

print("-" * 50)
print(f"True Escape Velocity: {v_true:.2f} m/s")
print(f"Yes, the calculation is possible. Smallest absolute error found is {abs_error:.2f} m/s.")

# The final answer format is specific
final_answer = f"Y{abs_error:.2f}"
# print(f"\n<<< {final_answer} >>>")