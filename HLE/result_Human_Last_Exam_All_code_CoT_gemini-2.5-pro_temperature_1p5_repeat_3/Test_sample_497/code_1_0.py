import sympy

def solve_ladder_capacitance():
    """
    This function solves for the terminating capacitance 'x' of a ladder network
    such that the equivalent capacitance is independent of the number of cells.
    The script follows the plan outlined above.
    """
    # --- Step 1: Define the symbolic variables ---
    # x represents the characteristic capacitance we are looking for.
    # c is the capacitance of the capacitors in each cell.
    x = sympy.Symbol('x', positive=True)
    c = sympy.Symbol('c', positive=True)

    print("Step 1: Formulate the equation for the characteristic capacitance.")
    print("The circuit is interpreted as a ladder of H-shaped cells, each with three capacitors of value 'c'.")
    print("Through circuit analysis (using symmetry), the input capacitance C_in for a cell terminated by a load C_load is found to be:")
    print("C_in = c * (c + C_load) / (3*c + 2*C_load)")
    print("\nFor the total capacitance to be independent of the number of cells, the input capacitance must equal the load capacitance.")
    print("This special value is the characteristic capacitance, which we will call 'x'.")
    print("So, we must solve the fixed-point equation:")
    print("x = c * (c + x) / (3*c + 2*x)")

    # --- Step 2: Rearrange into a quadratic equation ---
    print("\nStep 2: Rearrange the equation into the standard quadratic form (a*x\u00b2 + b*x + c = 0).")
    # x * (3*c + 2*x) = c * (c + x)
    # 3*c*x + 2*x**2 = c**2 + c*x
    # 2*x**2 + 2*c*x - c**2 = 0
    equation = sympy.Eq(2 * x**2 + 2 * c * x - c**2, 0)
    print("The resulting quadratic equation is:")
    print(f"2*x\u00b2 + (2*c)*x - c\u00b2 = 0")

    # --- Step 3: Solve the quadratic equation ---
    print("\nStep 3: Solve the quadratic equation for x.")
    solutions = sympy.solve(equation, x)
    print(f"The mathematical solutions for x are: {solutions[0]} and {solutions[1]}")

    # --- Step 4: Select the physically valid solution ---
    print("\nStep 4: Choose the physically valid solution.")
    # Since capacitance must be positive, we select the positive root.
    positive_solution = None
    for sol in solutions:
        # Check the sign by substituting a positive value like 1 for c
        if sol.subs(c, 1) > 0:
            positive_solution = sol
            break
            
    print("Capacitance must be a positive physical quantity, so we choose the positive root.")
    print(f"The required value for the capacitor x is: x = {positive_solution}")

    # --- Step 5: Output the final equation with its numerical components ---
    print("\nThe final equation for x can be constructed from the following numbers and operations:")
    # The expression is c*(-1 + sqrt(3))/2
    val_c = 'c'
    num_1 = 3
    num_2 = 1
    num_3 = 2
    print(f"x = {val_c} * (sqrt({num_1}) - {num_2}) / {num_3}")


if __name__ == '__main__':
    solve_ladder_capacitance()