import sympy

def solve_amplitude_equation():
    """
    This function solves the equation for the generating amplitude c1.
    The equation is derived from the analysis of the van der Pol oscillator,
    under the condition that the two amplitude components are equal (c1 = c2).
    """
    # Define the symbolic variable for the amplitude c1
    c1 = sympy.Symbol('c1')

    # The equation for the generating amplitude c1 is c1^2 - 2 = 0.
    # We can represent this with its coefficients.
    coeff_c1_squared = 1
    constant_term = -2
    
    # Form the equation
    equation = coeff_c1_squared * c1**2 + constant_term

    # Solve the equation for c1
    solutions = sympy.solve(equation, c1)

    # Find the first positive root
    positive_root = None
    for sol in solutions:
        if sol > 0:
            positive_root = sol
            break

    # Output the results, including the numbers in the final equation
    print("The equation for the generating amplitude c1 in the case c1 = c2 is:")
    # The format below explicitly shows each number in the equation
    print(f"{coeff_c1_squared} * c1**2 + ({constant_term}) = 0")
    
    if positive_root is not None:
        print(f"\nThe first positive root is c1 = {positive_root}")
        print(f"The numerical value is approximately: {positive_root.evalf()}")
    else:
        print("\nNo positive root was found.")

if __name__ == '__main__':
    solve_amplitude_equation()