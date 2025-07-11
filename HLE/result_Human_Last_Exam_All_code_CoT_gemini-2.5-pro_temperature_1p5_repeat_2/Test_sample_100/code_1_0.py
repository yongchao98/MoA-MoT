import sympy as sp

def solve_integral():
    """
    This function calculates the value of the integral by summing the four
    component series derived from the analytical solution.
    It prints the breakdown of the calculation as requested.
    """
    # Define symbols
    pi = sp.pi
    
    # Part 1: Sum from the p^7 term in the numerator
    # The resulting series is 7! * zeta(8)
    s1_val = sp.factorial(7) * sp.zeta(8)
    
    # Part 2: Sum from the p term in the numerator
    # The resulting series is zeta(2)
    s2_val = sp.zeta(2)
    
    # Part 3: Sum from the sinh(p/4) term in the numerator
    # The resulting series is 2 - pi/2
    s3_val = 2 - pi/2
    
    # Part 4: Sum from the p*exp(-p) term in the numerator
    # The resulting series is zeta(2) - 1
    s4_val = sp.zeta(2) - 1
    
    # Total value is the sum of these four parts
    total_integral = s1_val + s2_val + s3_val + s4_val
    
    # Simplify the final expression
    simplified_integral = sp.simplify(total_integral)
    
    print("The definite integral is computed by breaking it down into four infinite series:")
    
    # Output each component of the sum
    # Equation part 1: 7! * zeta(8)
    print(f"\nPart 1 = 7! * zeta(8)")
    print(f"       = {sp.factorial(7)} * ({pi**8}/9450)")
    print(f"       = {s1_val}")
    
    # Equation part 2: zeta(2)
    print(f"\nPart 2 = zeta(2)")
    print(f"       = {s2_val}")

    # Equation part 3: 2 - pi/2
    print(f"\nPart 3 = 2 - pi/2")
    print(f"       = {s3_val}")

    # Equation part 4: zeta(2) - 1
    print(f"\nPart 4 = zeta(2) - 1")
    print(f"       = {s4_val}")

    # Display the final combined equation
    print("\nThe final result is the sum of these four parts:")
    print(f"Result = ({s1_val}) + ({s2_val}) + ({s3_val}) + ({s4_val})")
    
    # Reconstruct the simplified expression string to show the "final equation"
    expr_terms = sp.Add.make_args(simplified_integral)
    term_list = [str(term) for term in expr_terms]
    final_eq_str = " + ".join(term_list).replace("+ -", "- ")
    
    print(f"Simplified Form: {final_eq_str}")

solve_integral()