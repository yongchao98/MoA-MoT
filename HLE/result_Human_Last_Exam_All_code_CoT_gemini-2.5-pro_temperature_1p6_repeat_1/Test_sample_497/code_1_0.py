import sympy

def solve_capacitance():
    """
    This function solves for the capacitance x such that the equivalent
    capacitance of the ladder network is independent of the number of cells.
    """
    # Step 1: Define the symbolic variables.
    # c is the capacitance of the capacitors in each cell.
    # C_load is the load capacitance terminating a cell.
    # C_in is the input capacitance of that terminated cell.
    c, C_load, C_in = sympy.symbols('c C_load C_in', positive=True)
    
    # Step 2: Define the recurrence relation for the input capacitance.
    # For a single H-shaped cell made of three capacitors 'c', loaded with C_load,
    # the input capacitance C_in can be found using circuit analysis (e.g., KCL).
    # The derived relation is:
    # C_in = c * (c + C_load) / (3*c + 2*C_load)
    # The condition for the total capacitance to be independent of the number of cells
    # is that the terminating capacitor x must be equal to the circuit's
    # characteristic capacitance, C_char. The characteristic capacitance is found
    # by setting the input capacitance equal to the load capacitance.
    # Let's call this value C_char.
    C_char = sympy.Symbol('C_char', positive=True)
    equation = sympy.Eq(C_char, c * (c + C_char) / (3*c + 2*C_char))
    
    print("Plan: Find the characteristic capacitance of the infinite ladder network.")
    print("Let the characteristic capacitance be C_char. The value of x should be C_char.")
    print("The characteristic capacitance is found by solving the fixed-point equation:\n")
    print(f"C_char = c * (c + C_char) / (3*c + 2*C_char)\n")

    # Step 3: Rearrange the equation into a quadratic form for C_char.
    # C_char * (3*c + 2*C_char) = c**2 + c*C_char
    # 3*c*C_char + 2*C_char**2 = c**2 + c*C_char
    # 2*C_char**2 + 2*c*C_char - c**2 = 0
    quadratic_eq = sympy.Eq(2*C_char**2 + 2*c*C_char - c**2, 0)
    print("This simplifies to the following quadratic equation for C_char:\n")
    print(f"{quadratic_eq}\n")
    
    # Step 4: Solve the quadratic equation for C_char.
    solutions = sympy.solve(quadratic_eq, C_char)
    
    # The solutions will contain a positive and a negative root. Capacitance must be positive.
    positive_solution = None
    for sol in solutions:
        # Check if the symbolic solution is positive given that c is positive.
        if sol.subs(c, 1) > 0:
            positive_solution = sol
            break
            
    print("Solving the quadratic equation gives two solutions. We choose the positive one as capacitance cannot be negative.\n")
    print(f"The positive solution for the characteristic capacitance is: C_char = {positive_solution}\n")
    
    # Step 5: The value of x is the characteristic capacitance.
    # The final equation highlights the numbers involved in the expression.
    # The expression is c * (sqrt(3) - 1) / 2
    number1 = 3
    number2 = 1
    number3 = 2
    
    print("--- Final Answer ---")
    print("The value of capacitor x must be equal to the characteristic capacitance.")
    print(f"The final equation for x in terms of c is:")
    print(f"x = c * (sqrt({number1}) - {number2}) / {number3}")

if __name__ == '__main__':
    solve_capacitance()