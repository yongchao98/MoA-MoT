import sympy

def solve_steady_state_pi0():
    """
    This function demonstrates the final step of the derivation for pi_0
    using the symbolic math library sympy.
    """
    # Define pi_0 and rho as symbolic variables
    pi_0 = sympy.Symbol('pi_0')
    rho = sympy.Symbol('rho')

    # The derivation shows that the normalization condition for the steady-state
    # probabilities simplifies to the equation: pi_0 * e^rho = 1.
    # We can represent this equation in sympy.
    # The term sympy.exp(rho) represents e^rho.
    # The number '1' is the right-hand side of the equation.
    final_equation = sympy.Eq(pi_0 * sympy.exp(rho), 1)
    
    # Solve the equation for the symbol pi_0
    solution = sympy.solve(final_equation, pi_0)
    
    # The solution is a list, so we take the first element
    pi_0_expression = solution[0]
    
    print("The final equation from the normalization condition is:")
    
    # We print each part of the equation: left side, equality, right side.
    # The number in the equation is 1.
    lhs = final_equation.lhs
    rhs = final_equation.rhs
    print(f"({lhs}) = {rhs}")
    
    print("\nSolving this equation for pi_0 yields the expression:")
    print(f"pi_0 = {pi_0_expression}")

if __name__ == "__main__":
    solve_steady_state_pi0()