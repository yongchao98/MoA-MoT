import sympy

def solve_ladder_capacitance():
    """
    Calculates the value of capacitor x such that the equivalent capacitance
    of the ladder network is independent of the number of cells.
    """
    # Define symbols for capacitances. 'c' and 'x' must be positive.
    c, x = sympy.symbols('c x', positive=True)
    C_char_half, C_load_half = sympy.symbols('C_char_half C_load_half', positive=True)

    print("Step 1: Model the circuit and set up the problem.")
    print("The circuit is a cascade of symmetric H-sections. The condition for capacitance to be independent of the number of cells (N) is that the capacitance of a finite ladder terminated by 'x' must equal the characteristic capacitance of an infinite ladder.")
    print("-" * 20)

    print("Step 2: Find the characteristic capacitance (C_char) using the half-circuit method.")
    print("The recurrence relation for the half-circuit's input capacitance (C_in_half) as a function of its load (C_load_half) is derived from its L-section structure (series 'c', shunt 'c/2'):")
    # C_in_half = (c^2/2 + c*C_load_half) / (3c/2 + C_load_half)
    recurrence_expr = (c**2/2 + c*C_load_half) / (sympy.Rational(3, 2)*c + C_load_half)
    print(f"C_in_half = {recurrence_expr}")
    
    # The characteristic capacitance is the fixed point of this relation (C_in_half = C_load_half = C_char_half)
    recurrence_eq = sympy.Eq(C_char_half, recurrence_expr.subs(C_load_half, C_char_half))

    # Solve for the characteristic capacitance of the half-ladder, C_char_half.
    # The quadratic equation 2*C^2 + c*C - c^2 = 0 yields two solutions.
    C_char_half_sols = sympy.solve(recurrence_eq, C_char_half)
    # We take the physically meaningful positive solution.
    C_char_half_val = [sol for sol in C_char_half_sols if sol.is_positive][0]

    # The characteristic capacitance of the full ladder is half of the half-circuit's.
    C_char_full = C_char_half_val / 2
    print(f"\nSolving the fixed-point equation gives the half-circuit characteristic capacitance: C_char_half = {C_char_half_val}")
    print(f"The characteristic capacitance of the full circuit is: C_char = C_char_half / 2 = {C_char_full}")
    print("-" * 20)
    
    print("Step 3: Calculate the capacitance of a 1-cell ladder terminated by x, C_1(x).")
    # For a 1-cell ladder terminated by x, the load on the half-circuit is 2*x.
    C_load_for_1_cell_half = 2*x
    C_1_half = recurrence_expr.subs(C_load_half, C_load_for_1_cell_half)
    
    # The total capacitance of the 1-cell ladder is C_1_half / 2
    C_1 = C_1_half / 2
    print(f"A 1-cell ladder's capacitance C_1(x) is found by setting C_load_half = 2x.")
    print(f"C_1(x) = {sympy.simplify(C_1)}")
    print("-" * 20)
    
    print("Step 4: Set C_1(x) equal to C_char and solve for x.")
    # For the capacitance to be independent of N, C_1(x) must equal C_char.
    final_eq = sympy.Eq(C_1, C_char_full)
    print(f"The equation to solve is: {sympy.simplify(C_1)} = {C_char_full}")
    
    # Solve the final equation for x.
    solution = sympy.solve(final_eq, x)
    final_x = solution[0]
    
    # Print the final result in the requested format.
    frac = final_x / c
    print("\n" + "="*40)
    print("The final result for x is:")
    print(f"x = {frac.p} * c / {frac.q}")
    print("Which can be written as:")
    print(f"x = {sympy.pretty(final_x, use_unicode=False)}")
    print("="*40)

if __name__ == '__main__':
    solve_ladder_capacitance()
