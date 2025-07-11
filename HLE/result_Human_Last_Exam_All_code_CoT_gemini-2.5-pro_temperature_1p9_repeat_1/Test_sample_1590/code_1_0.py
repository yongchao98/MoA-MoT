import sympy as sp

def solve_rod_sliding_angle():
    """
    This function provides the final expression for the angle at which the
    rod begins to slide, based on a physics derivation.
    """

    # Define the symbolic variables from the problem description.
    # L: total length of the rod
    # ell: distance from the pivot to the center of mass (denoted as l in the text)
    # mu: coefficient of static friction
    # theta: angle of the rod with the horizontal
    L, ell, mu = sp.symbols('L, ell, mu', positive=True, real=True)
    theta = sp.Symbol('theta')

    # After performing the physics derivation outlined in the plan,
    # the condition f = mu * N leads to the following expression for tan(theta):
    # tan(theta) = (mu * L**2) / (L**2 + 36 * ell**2)
    
    tan_theta_value = (mu * L**2) / (L**2 + 36 * ell**2)

    # The expression for the angle theta is the arctangent of this value.
    theta_expression = sp.atan(tan_theta_value)

    # Create the final equation object for clear printing.
    final_equation = sp.Eq(theta, theta_expression)

    # Print the final result in a human-readable format.
    print("The expression for the angle theta at which the rod begins to slide is:")
    # The sp.pretty() function formats the equation nicely.
    print(f"\n{sp.pretty(final_equation, use_unicode=False)}\n")

    # As requested, output the numbers present in the final derived expression.
    # The expression is arctan((mu * L**2) / (L**2 + 36 * ell**2)).
    # The constants are the numerical coefficients of the terms.
    print("The numerical constants in the expression are:")
    print("The coefficient of L**2 in the numerator is implicitly: 1")
    print("The coefficient of L**2 in the denominator is implicitly: 1")
    print("The coefficient of ell**2 in the denominator is: 36")

solve_rod_sliding_angle()