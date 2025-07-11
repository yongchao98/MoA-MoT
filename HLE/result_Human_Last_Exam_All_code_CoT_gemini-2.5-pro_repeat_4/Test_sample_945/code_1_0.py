import sympy as sp

def derive_and_print_critical_speed():
    """
    Derives and prints the formula for the critical speed of an oversteering vehicle
    using the linear single-track model and symbolic mathematics.
    """
    # Define symbolic variables for the vehicle parameters.
    # All are defined as positive real numbers.
    a, b, c_f, c_r, m, I, v = sp.symbols('a b c_f c_r m I v', positive=True, real=True)
    
    # --- Step 1: Explain the setup ---
    print("Derivation of Critical Speed for an Oversteering Vehicle")
    print("--------------------------------------------------------")
    print("The analysis is based on the linear single-track (bicycle) model.")
    print("The vehicle's lateral dynamics can be described by a state-space equation dx/dt = A * x.")
    print("The system becomes unstable when the determinant of the state matrix A is no longer positive.")
    print("The critical speed, v_crit, is the speed 'v' at which det(A) = 0.\n")

    # --- Step 2: Define the state matrix A ---
    print("The state matrix A for the linear single-track model is defined as:")
    # Elements of the state matrix A
    A11 = -(c_f + c_r) / (m * v)
    A12 = -1 - (a * c_f - b * c_r) / (m * v**2)
    A21 = (b * c_r - a * c_f) / I
    A22 = -(a**2 * c_f + b**2 * c_r) / (I * v)
    A = sp.Matrix([[A11, A12], [A21, A22]])
    sp.pprint(A)
    print("\n")
    
    # --- Step 3: Calculate the determinant of A ---
    print("Next, we calculate the determinant of A. For the system to be stable, det(A) must be > 0.")
    det_A = sp.det(A)
    
    # SymPy's det() result can be complex, so we simplify it.
    det_A_simplified = sp.simplify(det_A)
    print("det(A) =")
    sp.pprint(det_A_simplified)
    print("\n")

    # --- Step 4: Solve for the critical speed ---
    print("We find the critical speed by setting det(A) = 0 and solving for v.")
    print("Note: A real solution for v_crit only exists for an oversteering vehicle, where (a*c_f - b*c_r) > 0.\n")
    
    # Create the equation det(A) = 0
    equation_for_v = sp.Eq(det_A_simplified, 0)
    print("Equation to solve: det(A) = 0")
    # A cleaner representation of the equation is:
    # (c_f*c_r*(a+b)**2)/(m*I*v**2) - (a*c_f - b*c_r)/I = 0
    # Let's print this rearranged form to be clear:
    term1 = (c_f * c_r * (a+b)**2) / (m * I * v**2)
    term2 = (a * c_f - b * c_r) / I
    sp.pprint(sp.Eq(term1, term2))
    print("\n")


    # Solve for v^2
    solutions_v_sq = sp.solve(equation_for_v, v**2)
    v_crit_sq = solutions_v_sq[0] # There is only one solution for v^2

    print("Solving the equation for v^2, we get:")
    print("v_crit^2 =")
    sp.pprint(v_crit_sq)
    print("\n")

    # --- Step 5: Final Formula ---
    print("Taking the positive square root gives the final formula for the critical speed, v_crit:")
    print("v_crit = sqrt( Numerator / Denominator )\n")
    
    # To fulfill the request "output each number in the final equation", 
    # we will explicitly print the numerator and denominator of the expression under the square root.
    numerator = c_f * c_r * (a + b)**2
    denominator = m * (a * c_f - b * c_r)
    print(f"Numerator:   {numerator}")
    print(f"Denominator: {denominator}")
    
    # --- Step 6: Format final answer for submission ---
    final_formula_str = f"sqrt(({numerator})/({denominator}))"
    return final_formula_str

# --- Main Execution Block ---
if __name__ == '__main__':
    final_answer = derive_and_print_critical_speed()
    # The final output format required by the user prompt
    print(f"\n<<<v_crit = {final_answer}>>>")