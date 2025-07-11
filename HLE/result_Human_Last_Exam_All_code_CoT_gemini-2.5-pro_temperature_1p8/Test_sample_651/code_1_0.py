import sympy as sp

def solve_limit():
    """
    This function calculates the limit of M(theta) as theta -> 0.
    """
    # Define theta as a positive symbolic variable
    theta = sp.Symbol('theta', positive=True)

    # From the geometric analysis, the supremum of the angle alpha, M(theta),
    # is found to be equal to theta / 2.
    M_theta = theta / 2

    print(f"The expression for the supremum of the angle, M(theta), is: {M_theta}")

    # Calculate the limit of M(theta) as theta approaches 0.
    limit_value = sp.limit(M_theta, theta, 0)
    
    # Print the final equation and the numbers in it.
    print(f"The final equation is: limit(M(theta)) as theta -> 0 = {limit_value}")
    # The numbers in the final equation are the limit value on the RHS
    # and the value theta approaches on the LHS.
    print(f"The number theta approaches is: 0")
    print(f"The result of the limit is: {limit_value}")


solve_limit()