import math

def solve():
    """
    Calculates the approximate probability that a 2D simple random walk,
    starting from (3000, 4000) and conditioned to avoid the origin,
    never visits the origin's four neighbors.

    This probability is approximately 1 / (2 * |x|^2) for a starting point x far from the origin.
    """
    # Starting point coordinates
    x1 = 3000
    x2 = 4000

    # Calculate the squared distance from the origin
    squared_distance = x1**2 + x2**2

    # The constant in the asymptotic formula is 1/2
    C = 0.5

    # Calculate the approximate probability
    probability = C / squared_distance
    
    # We are asked for an approximate answer with two significant digits.
    # The final print out still has all the numbers used in the final equation.
    
    print(f"The starting point is ({x1}, {x2}).")
    print(f"The squared distance from the origin is |x|^2 = {x1}^2 + {x2}^2 = {squared_distance}.")
    print(f"The probability is approximated by the formula P(x) = C / |x|^2, where C = {C}.")
    print(f"P({x1},{x2}) = {C} / {squared_distance}")
    print(f"The calculated probability is: {probability:.2g}")

solve()
