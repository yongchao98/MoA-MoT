import math

# This script solves for the mixed-strategy equilibrium probability 'p' and computes the final requested value.

# Set the value of N for the game.
N = 8

def equation_to_solve(p, N):
    """
    Defines the root-finding equation f(p) = 0 derived from the equilibrium condition.
    The equation is: (1 - p/N)^(3N) - (1-p) = 0
    """
    # A check to prevent math domain errors, though p is expected to be in [0, 1].
    if p >= N or p < 0:
        return float('inf') 
    return (1.0 - p / N)**(3 * N) - (1.0 - p)

def solve_equilibrium_p(N):
    """
    Finds the non-zero root of the equilibrium equation using a bisection method.
    """
    # From analyzing the function, the non-zero root for N=8 lies between 0 and 1.
    # We can start with bounds [0.1, 0.99] to ensure we find the non-zero solution.
    low = 0.1
    high = 0.99
    
    # Iterate for 100 times to achieve high precision.
    for _ in range(100):
        mid = (low + high) / 2.0
        if equation_to_solve(mid, N) < 0:
            low = mid
        else:
            high = mid
    # The root is the converged value.
    return high

# Solve for p with N=8
p = solve_equilibrium_p(N)

# The problem asks for the value of floor(10000 * (1-p)).
# We calculate each part of this expression.
constant_multiplier = 10000
one_minus_p = 1.0 - p
product_value = constant_multiplier * one_minus_p
final_answer = math.floor(product_value)

# Print the numbers involved in the final calculation as requested.
print(f"The equation to be evaluated is: floor(A * (1-p))")
print(f"A = {constant_multiplier}")
print(f"p ≈ {p:.6f}")
print(f"1 - p ≈ {one_minus_p:.6f}")
print(f"The result is floor({constant_multiplier} * {one_minus_p:.6f}) = floor({product_value:.4f}) = {final_answer}")

<<<474>>>