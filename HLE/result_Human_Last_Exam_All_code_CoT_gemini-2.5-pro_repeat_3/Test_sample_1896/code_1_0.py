import numpy as np
from scipy.integrate import quad
import math

# Define the sinc function to avoid division by zero at x=0
def sinc(x):
    if x == 0:
        return 1.0
    else:
        return np.sin(x) / x

# Define the integrand for I_n, which is the product of sinc functions
def integrand(x, n):
    product = 1.0
    # The product is from k=1 to n
    for k in range(1, n + 1):
        product *= sinc(x / k)
    return product

# Value of pi/2 for comparison
pi_half = np.pi / 2

print("This script evaluates the Borwein integrals I_n for n from 1 to 9.")
print("It calculates the numerical value of the integral and its difference from pi/2.")
print("-" * 60)
print(f"Reference value of pi/2 = {pi_half:.15f}")
print("-" * 60)
print(f"{'n':<5}{'I_n (Numerical Value)':<25}{'Difference (I_n - pi/2)':<25}")
print("-" * 60)

# We will evaluate the integral for n from 1 up to 9 to see the behavior.
# The first deviation from pi/2 is known to occur at n=8.
for n in range(1, 10):
    # The quad function from SciPy performs numerical integration.
    # We integrate from 0 to infinity (np.inf).
    # The 'args' parameter is used to pass the current value of n to the integrand.
    # It returns the result and an estimate of the absolute error.
    result, error = quad(integrand, 0, np.inf, args=(n,))
    
    # Calculate the difference from the expected value for n <= 7
    difference = result - pi_half

    # Print the results in a formatted table
    print(f"{n:<5}{result:<25.15f}{difference:<25.15e}")

print("-" * 60)