import math

def solve_asymptotic_ode():
    """
    This function calculates the coefficient for the asymptotic solution of the ODE
    and prints the final analytical expression.
    """
    # From the analytical derivation, the coefficient A is determined by A^3 = -6.
    # So, A = - (6)^(1/3).
    A = -math.pow(6, 1.0/3.0)

    # Round the coefficient to two decimal places as requested.
    A_rounded = round(A, 2)

    # The problem asks to output each number in the final equation.
    # We will print the full expression for y(x).
    
    print("The derived analytical expression that approximates the solution for large x is:")
    print(f"y(x) = {A_rounded} / x")

solve_asymptotic_ode()