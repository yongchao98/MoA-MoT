import math

# This script calculates the number of iterations needed for the Schwarz
# Relaxation Method to converge for the 1D wave equation up to a time T.

# --- Problem Parameters ---
# T: The final time instant up to which convergence is required.
T = 2.5
# c: The propagation speed of the wave.
c = 1.0
# M: The size of the overlap region between subdomains (M = b - a).
M = 0.2

# --- Calculation ---
# The number of iterations 'K' is determined by the ratio of the total desired
# time 'T' to the time it takes for a wave to cross the overlap region 'M/c'.
# Since iterations are discrete, we take the ceiling of the result.

num_iterations = math.ceil(T * c / M)
product_T_c = T * c

# --- Output ---
# The final result is printed, showing the formula and the values used
# in the calculation as requested.
print("To find the number of iterations (K) for convergence up to time T, we use the formula:")
print("K = ceil(T * c / M)")
print("\nGiven the values:")
print(f"  Target time T = {T}")
print(f"  Propagation speed c = {c}")
print(f"  Overlap size M = {M}")
print("\nThe calculation is as follows:")
print(f"  K = ceil(({T} * {c}) / {M})")
print(f"  K = ceil({product_T_c} / {M})")
print(f"  K = ceil({product_T_c / M})")
print(f"  K = {int(num_iterations)}")