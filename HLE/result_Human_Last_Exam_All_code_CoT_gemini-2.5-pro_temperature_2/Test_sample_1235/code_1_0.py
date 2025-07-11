import sympy

def solve_amplitude_equation():
    """
    This function solves the equation for the generating amplitude of the
    van der Pol oscillator.
    """
    # Define the symbolic variable for the amplitude
    c = sympy.Symbol('c')

    # The equation for the non-trivial generating amplitude is c^2 = 4.
    # We define the numbers in the equation to print them explicitly as requested.
    power = 2
    value = 4
    
    # Formulate the equation in sympy
    equation = sympy.Eq(c**power, value)
    
    print("Assuming the problem describes the van der Pol oscillator, the equation for the non-trivial generating amplitude 'c' is:")
    
    # Output each number in the final equation
    print(f"c**{power} = {value}")

    # Solve the equation for real roots
    solutions = sympy.solve(equation, c, domain=sympy.S.Reals)

    # Find the first positive root c1
    c1 = None
    for s in solutions:
        if s > 0:
            c1 = s
            break
            
    if c1 is not None:
        print(f"\nThe first positive root of the equation is:")
        print(c1)
    else:
        print("\nNo positive root found.")

solve_amplitude_equation()
