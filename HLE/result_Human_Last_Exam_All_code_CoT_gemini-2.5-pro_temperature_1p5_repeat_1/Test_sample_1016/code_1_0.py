import math

# --- Parameters of the problem ---

# T: The final time instant up to which convergence is required (in seconds)
T = 10.0

# c: The propagation speed of the wave (in meters/second)
c = 2.0

# M: The size of the overlap between subdomains (in meters)
M = 1.0

# --- Calculation ---

# The convergence of the Schwarz method for the wave equation depends on how many
# times information can travel across the overlap region.
# Each iteration extends the convergence time by M/c.
# We need to find the smallest integer N such that N * (M/c) >= T.
# This leads to the formula: N = ceil(T * c / M).

# Calculate the value inside the ceiling function
value = T * c / M

# Calculate the number of iterations by taking the ceiling
# The result must be an integer.
num_iterations = math.ceil(value)

# --- Output the result ---

print("Problem parameters:")
print(f"  Final Time (T) = {T}")
print(f"  Wave Speed (c) = {c}")
print(f"  Overlap Size (M) = {M}")
print("\nCalculation of the number of iterations (N):")
print(f"  N = ceil(T * c / M)")
print(f"  N = ceil({T} * {c} / {M})")
print(f"  N = ceil({value})")
print(f"  N = {int(num_iterations)}")

# Final answer in the requested format
# print(f"\n<<<{int(num_iterations)}>>>")