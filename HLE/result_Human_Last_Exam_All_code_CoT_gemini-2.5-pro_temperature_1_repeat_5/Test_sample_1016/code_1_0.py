import math

# Description of the problem parameters
# T: The time instant up to which the method needs to converge.
# c: The propagation speed of the wave.
# M: The size of the overlap between the two subdomains (M = b - a).

# --- Example Parameters ---
# Let's consider a domain from x=0 to x=L=12.
# Wave propagation speed
c = 1.0
# Let the subdomains be Omega1 = [0, 8] and Omega2 = [4, 12].
# This means b=8 and a=4.
# The size of the overlap is M = b - a.
M = 8.0 - 4.0
# We want the solution to be converged up to the time instant T.
T = 20.0
# --------------------------

# The number of iterations 'k' for the Schwarz method to converge for the
# wave equation up to time T is determined by the time it takes for
# information to propagate across the overlap region in the iterative process.
#
# In each full iteration, information effectively travels across the overlap
# of size M and back, covering a "time" of 2*M/c.
# Therefore, after k iterations, the solution is converged up to time k * (2*M/c).
#
# We need k * (2*M/c) >= T.
# Solving for k: k >= T * c / (2 * M).
# Since k must be an integer, we take the ceiling.
# k = ceil(T * c / (2 * M))

# Perform the calculation
numerator = T * c
denominator = 2 * M
result_float = numerator / denominator
iterations = math.ceil(result_float)

# Output the results step-by-step
print("Calculation for the number of Schwarz iterations 'k'")
print("-" * 50)
print(f"Given Parameters:")
print(f"  Final Time (T)        = {T}")
print(f"  Wave Speed (c)        = {c}")
print(f"  Overlap Size (M)      = {M}")
print("-" * 50)
print("The formula for the number of iterations is: k = ceil(T * c / (2 * M))")
print("\nSubstituting the given values:")
print(f"k = ceil({T} * {c} / (2 * {M}))")
print(f"k = ceil({numerator} / {denominator})")
print(f"k = ceil({result_float:.4f})")
print(f"k = {iterations}")
print("-" * 50)
print(f"The method needs {iterations} iterations to converge up to time T = {T}.")
