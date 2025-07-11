import math

def solve_alpha():
    """
    Calculates the value of alpha based on the derived formula.
    """
    T = math.log(10)

    # The derived formula for alpha
    # alpha = (2 * (1 - exp(-3T))) / (3 * (1 - exp(-2T))) * exp(8T)

    c1 = 2
    c2 = -3
    c3 = 3
    c4 = -2
    c5 = 8

    term1 = 1 - math.exp(c2 * T)
    term2 = 1 - math.exp(c4 * T)
    term3 = math.exp(c5 * T)

    alpha = (c1 * term1) / (c3 * term2) * term3

    print(f"The final equation for alpha is: alpha = ({c1} * (1 - exp({c2}*T))) / ({c3} * (1 - exp({c4}*T))) * exp({c5}*T)")
    print(f"With T = ln(10), the calculated value for alpha is:")
    print(alpha)

solve_alpha()