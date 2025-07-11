import sympy as sp

def calculate_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) based on the problem statement.
    It uses the derived system of ODEs for the coordinates in the given frame and calculates
    the rate of change of the angle in their polar representation.
    """
    # Define symbols for the constant c, the Gaussian curvature K, and the angle theta.
    # These represent c, K(gamma(t)), and theta(t) respectively.
    c = sp.Symbol('c')
    K = sp.Symbol('K')
    theta = sp.Symbol('theta')

    # The coordinates (A, B) in the given frame can be written in polar form.
    # Since the radius r(t) cancels out in the expression for theta'(t), we can set r=1.
    A = sp.cos(theta)
    B = sp.sin(theta)
    
    # From the derivation, the system of ODEs for (A, B) is:
    # A'(t) = -(K/c) * B(t)
    # B'(t) = c * A(t)
    # We define the expressions for the derivatives Adot and Bdot.
    Adot = -(K / c) * B
    Bdot = c * A

    # The formula for the derivative of the angle is theta' = (A*B' - B*A') / (A^2 + B^2).
    # We substitute our expressions and simplify.
    theta_prime = sp.simplify((A * Bdot - B * Adot) / (A**2 + B**2))

    # The problem asks to output the final equation. We will print the
    # simplified result in a standard mathematical notation.
    print("The final expression for theta'(t) is:")
    sp.pretty_print(theta_prime)

    # To fulfill the request to output each component of the final equation,
    # we can break down the expression.
    term1 = c * sp.cos(theta)**2
    term2 = (K / c) * sp.sin(theta)**2
    
    print("\nThe equation is composed of two terms added together:")
    print("Term 1:", term1)
    print("Term 2:", term2)

calculate_theta_prime()