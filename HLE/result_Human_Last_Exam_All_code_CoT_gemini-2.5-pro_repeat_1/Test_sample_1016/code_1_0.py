import math

# Define the parameters for the wave equation and domain decomposition.
# These are example values.
# c: wave propagation speed (e.g., in meters/second)
c = 343.0
# T: The time instant up to which convergence is required (in seconds)
T = 0.05
# a: Start of the right subdomain (in meters)
a = 4.0
# b: End of the left subdomain (in meters)
b = 5.0

# --- Calculation ---

# 1. Calculate the size of the overlap region M.
M = b - a

# 2. Derive the condition for convergence.
# The method converges when the information has propagated a distance of c*T.
# Each iteration allows information to cross the overlap region of size M.
# So, for K iterations, the condition is: K * M >= c * T
# This means K must be at least (c * T) / M.

# 3. Calculate the minimum number of integer iterations.
# Since the number of iterations K must be an integer, we take the
# ceiling of the division.
physical_distance = c * T
num_iterations_float = physical_distance / M
num_iterations = math.ceil(num_iterations_float)

# --- Output the results ---

print("Problem: Find the number of Schwarz iterations for the 1D wave equation.")
print("-" * 60)
print(f"Given Parameters:")
print(f"  Wave speed (c)              : {c}")
print(f"  Convergence time (T)        : {T}")
print(f"  Left subdomain end (b)      : {b}")
print(f"  Right subdomain start (a)   : {a}")
print("")

print("Step 1: Calculate the overlap size M = b - a")
print(f"  M = {b} - {a}")
print(f"  M = {M}")
print("")

print("Step 2: Calculate the required number of iterations K.")
print("  The formula is K = ceil((c * T) / M)")
print(f"  K = ceil(({c} * {T}) / {M})")
print(f"  K = ceil({physical_distance} / {M})")
print(f"  K = ceil({num_iterations_float})")
print(f"  K = {num_iterations}")
print("-" * 60)

print(f"\nThe method requires {num_iterations} iterations to converge up to time T = {T}s.")
