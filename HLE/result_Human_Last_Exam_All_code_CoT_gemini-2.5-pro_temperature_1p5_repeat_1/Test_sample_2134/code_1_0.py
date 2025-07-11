import math

def solve_trajectory():
    """
    Calculates the final position of the particle based on the derived trajectory equation.
    The final position x(t) at t=2*sqrt(3) is given by the expression:
    x(2*sqrt(3)) = (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
    """

    # Value of sqrt(3)
    sqrt3 = math.sqrt(3)

    # Calculate the values inside the cubic roots
    val1 = 18 - 6 * sqrt3
    val2 = 18 + 6 * sqrt3

    # Calculate each term of the sum
    term1 = val1**(1/3)
    term2 = val2**(1/3)

    # Calculate the final result
    result = term1 + term2
    
    # Print the equation and the final calculated value
    print(f"The equation for the final position is: x = (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)")
    print(f"The first term evaluates to: ({val1})^(1/3) = {term1}")
    print(f"The second term evaluates to: ({val2})^(1/3) = {term2}")
    print(f"The final value of x(2*sqrt(3)) is the sum of these terms: {term1} + {term2} = {result}")

solve_trajectory()