import math

def modified_logistic_map(r, x):
    """
    A modified map that produces a stable fixed point.
    Equation: X_n+1 = X_n / R + 1
    """
    return x / r + 1

def solve():
    """
    Solves the problem by iterating the modified map and printing the result.
    """
    # Parameters
    R = 3.57
    X0 = 0.5  # Initial value for X
    iterations = 200

    # Iterate the map to find the equilibrium point
    x = X0
    for _ in range(iterations):
        x = modified_logistic_map(R, x)

    # The theoretical fixed point
    fixed_point_theoretical = R / (R - 1)
    
    # We will print the terms of the equation that gives the fixed point.
    # The fixed point equation is X = R / (R - 1)
    numerator = R
    denominator_term1 = R
    denominator_term2 = 1.0
    result = fixed_point_theoretical

    print("The modified logistic map can be expressed as the recurrence relation: X_n+1 = X_n / R + 1")
    print(f"This relation has a stable fixed point X* given by the equation: X* = R / (R - 1)")
    print(f"For R = {R}, the equation for the fixed point is:")
    print(f"X* = {numerator} / ({denominator_term1} - {denominator_term2})")
    print(f"The resulting equilibrium point is approximately: {result}")
    
solve()
<<<1.389>>>