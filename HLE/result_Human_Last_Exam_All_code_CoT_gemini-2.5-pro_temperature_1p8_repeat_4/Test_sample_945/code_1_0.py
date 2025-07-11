import sympy as sp

def derive_and_display_critical_speed():
    """
    This function symbolically derives and displays the formula for the critical speed
    of an oversteering vehicle using the linear single-track model.
    It explains the process and prints the final equation in a clear format.
    """
    
    # --- Introduction ---
    print("Derivation of Critical Speed for an Oversteering Vehicle")
    print("---------------------------------------------------------")
    print("This script derives the critical speed using stability analysis of the linear single-track model.")
    print("The parameters used are:")
    print("  a: distance from CG to front axle")
    print("  b: distance from CG to rear axle")
    print("  c_f: cornering stiffness of the front axle")
    print("  c_r: cornering stiffness of the rear axle")
    print("  m: vehicle mass")
    print("  I: vehicle moment of inertia")
    print("  v: vehicle forward speed\n")

    # --- Symbolic Derivation ---
    print("--- Derivation Steps ---")
    # Step 1: Define symbolic variables
    a, b, cf, cr, m, I, v = sp.symbols('a b c_f c_r m I v', positive=True)
    v_crit = sp.Symbol('v_crit')

    # The stability of the system is determined by the determinant of its state matrix A.
    # The system becomes unstable when the determinant is no longer positive.
    # The determinant is given by:
    # Det(A) = (c_f*c_r*(a+b)**2) / (m*I*v**2) + (b*c_r - a*c_f) / I
    
    print("1. The stability of the vehicle's lateral motion is governed by the sign of the determinant of its state matrix 'A'.")
    print("   Det(A) = (c_f*c_r*(a+b)**2)/(m*I*v**2) - (a*c_f - b*c_r)/I")
    
    print("\n2. An oversteering vehicle is defined by the condition: a*c_f > b*c_r.")
    print("   This makes the term 'a*c_f - b*c_r' positive.")

    print("\n3. The critical speed is the speed 'v' at which the vehicle becomes unstable, which occurs when Det(A) = 0.")
    print("   Setting Det(A) = 0 and solving for v gives:\n")
    print("   (c_f*c_r*(a+b)**2) / (m*I*v_crit**2) = (a*c_f - b*c_r) / I")
    
    # Step 4: Solve for the critical speed
    # Rearranging the equation to solve for v_crit^2
    # v_crit^2 = (c_f * c_r * (a+b)^2 * I) / (m * I * (a*c_f - b*c_r))
    
    numerator_expr = cf * cr * (a + b)**2
    denominator_expr = m * (a * cf - b * cr)
    critical_speed_expr = sp.sqrt(numerator_expr / denominator_expr)
    
    final_equation = sp.Eq(v_crit, critical_speed_expr)

    # --- Final Result ---
    print("\n--- Final Equation for Critical Speed ---")
    print("The resulting equation shows each symbolic parameter's place in the formula:")
    
    # Print the "raw" formula as a string, showing each term
    print(f"\nv_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r)) )")
    
    # Use sympy's pretty print for a nicely formatted symbolic representation
    print("\nSymbolic representation of the final equation:")
    sp.pprint(final_equation, use_unicode=True)
    
    # Return the string representation for the final answer block
    return str(critical_speed_expr)

if __name__ == '__main__':
    # Execute the derivation and printing function
    final_answer_str = derive_and_display_critical_speed()
    # Present the final answer in the required format
    print(f"\n<<<{final_answer_str}>>>")
