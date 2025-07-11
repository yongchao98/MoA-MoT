import math

def solve_control_problem():
    """
    Calculates the value of the given expression based on the control problem.

    The problem is a linear control system:
    x_1'(t) = x_1(t) + 2x_2(t) + 3u_1
    x_2'(t) = 2x_1(t) + 4u_2
    with u_1 = (2/3)u_2.

    The boundary condition at t = 1/2 is:
    (2/3)x_1(1/2) - (1/3)x_2(1/2) = 4

    The eigenvalues are given as:
    lambda_1 = (1 + sqrt(17))/2
    lambda_2 = (1 - sqrt(17))/2

    The expression to find is:
    ( (2/3)lambda_1*e^(lambda_1/2) - (1/3)e^(lambda_1/2) ) * x_2(0) 
    - (2/3)lambda_2*e^(lambda_2/2) - (10/3)e^(lambda_1/2)

    A key step is to make a reasonable assumption to make the problem solvable.
    In many controllability problems, it is assumed that the system starts from
    the origin, i.e., x(0) = 0. This implies x_1(0) = 0 and x_2(0) = 0.
    Under this assumption, the first term of the expression becomes zero.
    """

    lambda_1 = (1 + math.sqrt(17)) / 2
    lambda_2 = (1 - math.sqrt(17)) / 2

    # Assumption: The initial state is the origin.
    x2_0 = 0

    # Calculate the components of the expression.
    # The expression is E = C1 * x2(0) + C2 + C3
    
    # Term 1: Coefficient of x2(0)
    c1 = (2/3 * lambda_1 - 1/3) * math.exp(lambda_1 / 2)
    term1 = c1 * x2_0
    
    # Term 2
    c2 = - (2/3) * lambda_2 * math.exp(lambda_2 / 2)
    
    # Term 3
    c3 = - (10/3) * math.exp(lambda_1 / 2)

    # Final value
    final_value = term1 + c2 + c3
    
    # The final equation is Value = term1 + term2 + term3
    print(f"The final equation is: ({c1}) * {x2_0} + ({c2}) + ({c3})")
    print(f"Breaking it down into constituent numerical parts:")
    print(f"First number (coefficient of x2(0)): {c1}")
    print(f"Value of x2(0) (assumed): {x2_0}")
    print(f"Second number: {c2}")
    print(f"Third number: {c3}")
    
    # Let's print the result of each main term.
    print(f"The equation evaluates to: {term1} + {c2} + {c3} = {final_value}")
    
    print(f"\nThe final calculated value is:")
    print(final_value)

solve_control_problem()