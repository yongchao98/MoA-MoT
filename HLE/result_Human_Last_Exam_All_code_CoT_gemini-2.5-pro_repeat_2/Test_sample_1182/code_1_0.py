def solve_curve_reduction():
    """
    This script explains the steps to find the number of double points
    in the stable reduction of the curve y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5.
    The actual calculations are mathematical and are explained in the comments.
    """

    # The original curve is C: y^2 = 8*x^5 + 4*x^4 + 4*x^3 + 1*x^2 + 8*x.
    # We find a new model C' over the 2-adic integers by the transformation
    # x = X/2, y = Y/2. This leads to the equation:
    # Y^2 = 1*X^5 + 1*X^4 + 2*X^3 + 1*X^2 + 16*X.
    
    # The stable reduction is the reduction of this model modulo 2.
    # Y^2 = X^5 + X^4 + X^2 (mod 2)
    # This can be factored as:
    # Y^2 = X^2 * (X^3 + X^2 + 1)
    
    # We print the coefficients of the factored form of the stable reduction equation.
    # The equation is 1*Y^2 = 1*X^2 * (1*X^3 + 1*X^2 + 0*X + 1).
    print("The equation for the stable reduction of the curve is Y^2 = X^2 * (X^3 + X^2 + 1).")
    print("The coefficients of the terms Y^2, X^2, X^3, X^2, 1 are respectively:")
    print(1) # for Y^2
    print(1) # for X^2
    print(1) # for X^3
    print(1) # for X^2 inside the parenthesis
    print(1) # for the constant term 1
    
    # As explained in the thought process, this reduced curve is a stable curve.
    # Its only singularity is an ordinary double point (a node) at (0,0).
    num_double_points = 1
    
    print("\nThe number of double points in the stable reduction is:")
    print(num_double_points)

solve_curve_reduction()