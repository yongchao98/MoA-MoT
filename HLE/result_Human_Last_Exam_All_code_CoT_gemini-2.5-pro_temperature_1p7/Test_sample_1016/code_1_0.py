import math

# --- User-defined parameters ---
# Please replace these example values with your specific problem parameters.
# c: Propagation speed of the wave
# T: The time instant up to which convergence is required
# M: The size of the overlap between subdomains (M = b - a)

c = 343.0  # Example: Speed of sound in air (m/s)
T = 0.1    # Example: 0.1 seconds
M = 1.0    # Example: 1.0 meter overlap

# --- Calculation ---
# The Schwarz Relaxation Method for the wave equation with absorbing boundary
# conditions at the interfaces converges in a finite number of iterations for any
# finite time interval [0, T].
#
# The convergence is determined by the time it takes for information to make a
# round trip across the overlap region of size M.
# Time for a round trip = 2 * M / c
#
# For the solution to be correct up to time T, the number of iterations N must satisfy:
# N * (Time for a round trip) >= T
# N * (2 * M / c) >= T
#
# Solving for N gives the formula:
# N >= (c * T) / (2 * M)
#
# Since N must be an integer, we take the ceiling of the result.

# Calculate the ratio (c*T) / (2*M)
value = (c * T) / (2 * M)

# The number of iterations is the ceiling of this value.
iterations = math.ceil(value)

# --- Output ---
print("This script calculates the number of iterations required for the Schwarz method to converge for the 1D wave equation.")
print("="*80)
print("Given parameters:")
print(f"  Propagation Speed (c) = {c}")
print(f"  Overlap Size (M)      = {M}")
print(f"  Final Time (T)        = {T}")
print("-"*80)
print("The formula for the required number of iterations (N) is:")
print("  N = ceil((c * T) / (2 * M))")
print("\nCalculation steps:")
print(f"  N = ceil(({c} * {T}) / (2 * {M}))")
print(f"  N = ceil({c * T} / {2 * M})")
print(f"  N = ceil({value})")
print("="*80)
print(f"Result: The method needs to perform {int(iterations)} iterations to converge up to time T={T}.")