import sympy

def find_fixed_point_coupling():
    """
    Calculates and displays the leading order expression for the fixed point
    coupling u* in phi^4 theory near d=4 dimensions.
    """
    # Define the symbols for the coupling and epsilon
    u = sympy.Symbol('u^*')
    epsilon = sympy.Symbol('epsilon')
    pi = sympy.pi

    print("The fixed point coupling u^* in phi^4 theory is determined by the condition that the beta function is zero: beta(u^*) = 0.")
    print(f"In d = 4 - epsilon dimensions, the one-loop beta function for the coupling u in L_int = (u/4!)*phi^4 is:")
    print(f"beta(u) = -epsilon * u + (3 / (16 * pi^2)) * u^2\n")

    # Define the coefficients from the beta function equation beta(u*) = 0
    # Equation: -epsilon*u + (c_num / (c_den * pi**2)) * u**2 = 0
    c_num = 3
    c_den = 16
    
    # Define the equation using sympy
    equation = -epsilon * u + (c_num / (c_den * pi**2)) * u**2
    
    print("Setting beta(u^*) = 0 gives the equation:")
    print(f"   - {epsilon} * {u} + ({c_num} / ({c_den} * pi^2)) * {u}^2 = 0")
    print("\nWe can factor out u^*:")
    print(f"   {u} * (-{epsilon} + ({c_num} / ({c_den} * pi^2)) * {u}) = 0")

    # Solve the equation for u
    solutions = sympy.solve(equation, u)

    # The solutions will be a list containing 0 and the non-trivial solution.
    trivial_solution = solutions[0]
    non_trivial_solution = solutions[1]
    
    print("\nThis equation has two fixed point solutions:")
    print(f"1. The trivial (Gaussian) fixed point: {u} = {trivial_solution}")
    print(f"2. The non-trivial (Wilson-Fisher) fixed point, obtained by solving:")
    print(f"   -{epsilon} + ({c_num} / ({c_den} * pi^2)) * {u} = 0")
    
    print("\nSolving for u^* gives:")
    print(f"   ({c_num} / ({c_den} * pi^2)) * {u} = {epsilon}")
    print(f"   {u} = {epsilon} * ({c_den} * pi^2) / {c_num}")

    print("\nTherefore, the leading order expression for the Wilson-Fisher fixed point coupling is:")
    # Pretty print the final symbolic expression
    final_expression = sympy.pretty(non_trivial_solution, use_unicode=False)
    # The final expression from sympy might be on multiple lines. We'll format it a bit.
    final_expression_single_line = str(non_trivial_solution).replace('**', '^')
    print(f"   {u} = {final_expression_single_line}")

if __name__ == '__main__':
    find_fixed_point_coupling()