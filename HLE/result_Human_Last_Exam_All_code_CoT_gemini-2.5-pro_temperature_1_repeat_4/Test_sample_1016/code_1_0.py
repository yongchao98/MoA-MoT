import math

# --- Problem Parameters ---
# These values are examples. You can change them to fit a specific scenario.

# c: The propagation speed of the wave (e.g., in m/s)
c = 343.0

# T: The final time up to which convergence is required (e.g., in seconds)
T = 0.1

# M: The size of the overlap region, M = b - a (e.g., in meters)
M = 10.0

# --- Calculation ---

# The Schwarz method for the wave equation converges for a time t after k iterations
# where t < k * (M / c).
# To have convergence up to time T, we need k * (M / c) >= T.
# Solving for k gives: k >= (c * T) / M

# Calculate the value of (c * T) / M
min_k_float = (c * T) / M

# Since the number of iterations k must be an integer, we take the ceiling
# of the result. The iteration count starts at 0, so k represents the number
# of full iteration cycles needed.
k = math.ceil(min_k_float)


# --- Output ---

print("Calculating the number of iterations for Schwarz method convergence.")
print("-" * 60)
print(f"Wave propagation speed (c): {c}")
print(f"Final simulation time (T):   {T}")
print(f"Overlap size (M):            {M}")
print("-" * 60)

print("The formula for the required number of iterations (k) is:")
print("k >= (c * T) / M\n")

print("Substituting the values into the equation:")
# The print statement below fulfills the requirement to output each number in the final equation
print(f"k >= ({c} * {T}) / {M}")
print(f"k >= {c * T} / {M}")
print(f"k >= {min_k_float:.4f}\n")

print("Since k must be an integer, we take the ceiling of this value:")
print(f"k = ceil({min_k_float:.4f})")
print(f"k = {k}\n")

print(f"Result: The method requires {k} iterations to converge up to time T = {T}.")
