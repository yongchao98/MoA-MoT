import sympy

def solve_ladder_capacitance():
    """
    Solves for the terminating capacitance x such that the equivalent capacitance
    of the ladder circuit is independent of the number of cells N.

    The final result for x is printed, showing each component of the equation.
    """
    # --- Step 1: Define symbolic variables ---
    # c: capacitance of the capacitors in a cell
    # C_k: equivalent capacitance of a ladder with k cells to the right
    # x: the terminating capacitance
    # We will use 'C_k' to represent the capacitance of the rest of the ladder
    # and 'C_k_minus_1' to represent the new total capacitance after adding one cell.
    c, C_k, C_k_minus_1 = sympy.symbols('c C_k C_{k-1}', positive=True, real=True)
    x = sympy.Symbol('x', positive=True, real=True)

    print("Step 1: Establishing the recurrence relation for the equivalent capacitance.")
    print("Let C_k be the equivalent capacitance of the circuit to the right of a given point.")
    
    # The load for a new cell is the parallel combination of the shunt capacitor 'c'
    # and the rest of the ladder 'C_k'.
    C_load = c + C_k
    
    # The new equivalent capacitance C_k_minus_1 is derived from the T-network of
    # two series 'c' capacitors and the shunt load C_load. The relation is:
    recurrence_formula = c * C_load / (c + 2 * C_load)
    
    # Substitute C_load to get the final recurrence relation
    final_recurrence = sympy.simplify(recurrence_formula.subs(C_load, c + C_k))
    print(f"The recurrence relation is C_{{k-1}} = {final_recurrence}")
    print("-" * 30)

    # --- Step 2: Find the characteristic capacitance (C_inf) ---
    print("Step 2: Finding the characteristic capacitance of an infinite ladder.")
    print("This is the 'fixed point' where adding a cell doesn't change the capacitance.")
    
    C_inf = sympy.Symbol('C_inf', positive=True, real=True)
    # We find the fixed point by setting C_k = C_{k-1} = C_inf
    equation = sympy.Eq(C_inf, final_recurrence.subs(C_k, C_inf))
    print("We solve the fixed-point equation:")
    sympy.pretty_print(equation)
    
    # Solve the quadratic equation for C_inf
    solutions = sympy.solve(equation, C_inf)
    
    # The characteristic capacitance must be a positive value.
    characteristic_capacitance = solutions[0] # Sympy takes care of the positive assumption
    
    print("\nThe positive solution for the characteristic capacitance is:")
    sympy.pretty_print(sympy.Eq(C_inf, characteristic_capacitance))
    print("-" * 30)

    # --- Step 3: Determine the value of x ---
    print("Step 3: Determining the value of the terminating capacitor x.")
    print("For the total equivalent capacitance to be independent of the number of cells,")
    print("the ladder must be terminated in its own characteristic capacitance.")
    print("This is achieved by setting the terminator x equal to C_inf.")
    
    x_value = characteristic_capacitance
    
    # Extracting the components of the answer for clear printing as requested.
    # The result is of the form c * (sqrt(A) - B) / D
    # For our result c*(-1 + sqrt(3))/2: A=3, B=1, D=2
    A = 3
    B = 1
    D = 2
    
    print("\n--- Final Answer ---")
    print("The value of the capacitor x is given by the equation:")
    print(f"x = c * (sqrt({A}) - {B}) / {D}")
    print("\nWhich simplifies to:")
    final_eq = sympy.Eq(x, x_value)
    sympy.pretty_print(final_eq)

solve_ladder_capacitance()