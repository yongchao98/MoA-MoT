import math

def get_true_force():
    """Calculates the gravitational force using standard floating-point arithmetic."""
    G = 6.67430e-11
    m_probe = 30.0  # kg

    # Pandora's properties in SI units
    R_core = 50 * 1000  # 50 km in meters
    R_outer = 1000 * 1000  # 1000 km in meters
    rho_core = 1.2 * 1000  # 1.2 metric tons/m^3 in kg/m^3
    rho_shell = 0.3 * 1000  # 0.3 metric tons/m^3 in kg/m^3
    
    # Probe's distance from center
    altitude = 500  # meters
    r = R_outer + altitude

    # Mass calculation
    M_core_mass = rho_core * (4/3) * math.pi * (R_core**3)
    V_shell = (4/3) * math.pi * (R_outer**3) - (4/3) * math.pi * (R_core**3)
    M_shell_mass = rho_shell * V_shell
    M_total = M_core_mass + M_shell_mass

    # Force calculation
    force = (G * M_total * m_probe) / (r**2)
    return force

def titan_multiply(frac1, frac2):
    """
    Multiplies two fractions according to Titan's rules.
    It performs cross-simplification before multiplication to keep intermediate values small.
    """
    n1, d1 = frac1
    n2, d2 = frac2

    # Check if input fractions are valid
    if not all(0 <= x <= 31 for x in [n1, d1, n2, d2]):
        raise ValueError("Input numerators/denominators must be between 0 and 31.")

    # Cross-simplify before multiplying
    common1 = math.gcd(n1, d2)
    n1_s, d2_s = n1 // common1, d2 // common1

    common2 = math.gcd(n2, d1)
    n2_s, d1_s = n2 // common2, d1 // common2
    
    # Multiply the simplified parts
    final_n = n1_s * n2_s
    final_d = d1_s * d2_s

    # Final check for overflow
    if final_n > 31 or final_d > 31:
        # This indicates the operation is invalid without re-approximation
        # For this problem, we've chosen fractions that avoid this
        raise ValueError(f"Resulting fraction {final_n}/{final_d} exceeds 5-bit limit.")
        
    return final_n, final_d

def solve_landing_problem():
    """
    Solves the Pandora landing problem using Titan's computational rules.
    """
    # 1. Define the components of the simplified formula F ≈ (4/3) * pi * C
    term1 = (4, 3)  # Represents 4/3

    # 2. Approximate C = G*m*rho_shell*R_outer ≈ 0.6007 with 3/5
    # C_val = 0.6007, C_approx = 3/5 = 0.6
    term_C = (3, 5) 

    # 3. Approximate pi ≈ 3.14159 with 25/8
    # This is chosen specifically because 25 simplifies with the 5 in C's denominator
    # and 8 simplifies with the 4 in the first term.
    # pi_val = 3.14159, pi_approx = 25/8 = 3.125
    term_pi = (25, 8)

    # 4. Perform the calculation step-by-step, reordering for simplification
    # Expression: (4/3) * (25/8) * (3/5)
    # Reorder to: (4/3) * (3/5) * (25/8)
    
    # Step 1: (4/3) * (3/5)
    intermediate_result = titan_multiply(term1, term_C) # (4,5)

    # Step 2: (4/5) * (25/8)
    final_result = titan_multiply(intermediate_result, term_pi) # (5,2)

    # 5. Calculate the error
    force_titan = final_result[0] / final_result[1]
    force_true = get_true_force()
    absolute_error = abs(force_titan - force_true)
    
    # 6. Print the final answer
    # The problem asks to output the equation with each number.
    print("Final Titan Equation:")
    # We show the reordered calculation that works
    print(f"{term1[0]} / {term1[1]} * {term_C[0]} / {term_C[1]} * {term_pi[0]} / {term_pi[1]} = {final_result[0]} / {final_result[1]}")
    
    print(f"\nCalculated Force (Titan): {force_titan:.3f} N")
    print(f"Theoretical Force: {force_true:.3f} N")
    print(f"Absolute Error: {absolute_error:.3f}")

    # Output the final answer in the required format
    print("\n" + "<<<" + f"Y{absolute_error:.3f}" + ">>>")

solve_landing_problem()