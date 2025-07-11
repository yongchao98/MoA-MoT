import math
from scipy.optimize import brentq

# N = 8 races
N = 8

# Define the polynomial function f(x) = 191x^24 - 192x^23 + 1
def f(x):
    """
    Polynomial equation derived from the equilibrium condition.
    """
    # Using math.pow for precision, though ** is usually fine
    return 191.0 * math.pow(x, 24) - 192.0 * math.pow(x, 23) + 1.0

# The probability p must be in (0, 1).
# x = 1 - p/N = 1 - p/8
# So x must be in (1-1/8, 1) = (7/8, 1) = (0.875, 1)
a = 0.875
b = 1.0

# Find the root of f(x) in the interval (a, b)
# brentq is a reliable root-finding algorithm
try:
    x_star = brentq(f, a, b)
    
    # Calculate p from the root x*
    p = N * (1 - x_star)
    
    # The quantity to calculate is floor(10000 * (1-p))
    quantity = 1 - p
    result = math.floor(10000 * quantity)
    
    # Print the intermediate values as requested
    print(f"For N = {N}, we solve the equation: 191*x^24 - 192*x^23 + 1 = 0")
    print(f"The non-trivial root is found in the interval ({a}, {b}).")
    print(f"Numerical solver found root x = {x_star:.8f}")
    
    f_at_x_star = f(x_star)
    print(f"Checking the function at the root: f(x) = {f_at_x_star:.2e}")

    print(f"This gives p = {N} * (1 - x) = {N} * (1 - {x_star:.8f}) = {p:.8f}")
    print(f"The problem asks for the value of floor(10000 * (1-p))")
    print(f"1 - p = {quantity:.8f}")
    print(f"10000 * (1 - p) = {10000 * quantity:.4f}")
    print(f"The final answer is floor({10000 * quantity:.4f}) = {result}")

except ValueError:
    print("Could not find a root in the given interval.")
    result = "Error"

print(f"<<<{result}>>>")
