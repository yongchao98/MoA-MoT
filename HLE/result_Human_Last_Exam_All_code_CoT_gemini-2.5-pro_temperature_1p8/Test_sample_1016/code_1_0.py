import math

# --- Parameters ---
# You can change these values to match your specific problem.
# c: propagation speed of the wave
c = 343.0  # m/s (e.g., speed of sound in air)
# T: The final time up to which the solution needs to be correct
T = 0.1  # seconds
# M: The size of the overlap region (M = b - a)
M = 5.0  # meters

# --- Calculation ---
# The number of iterations 'k' required for the Schwarz method to converge
# for the 1D wave equation up to time T is determined by the time it takes
# for information to make a round trip across the overlap region.
#
# Time for a one-way trip across the overlap = M / c
# Time for a round trip (information exchange) = 2 * M / c
#
# Each iteration completes one cycle of information exchange, so it advances
# the convergence time by 2 * M / c. To be converged up to time T, we need:
# k * (2 * M / c) >= T
#
# Solving for k gives:
# k >= (c * T) / (2 * M)
#
# Since the number of iterations must be an integer, we take the ceiling.

print("This script calculates the number of iterations for Schwarz method convergence.")
print("-" * 70)

if c <= 0 or T < 0 or M <= 0:
    print("Error: Parameters c, T, and M must be positive.")
else:
    # Calculate the expression (c * T) / (2 * M)
    numerator = c * T
    denominator = 2 * M
    result = numerator / denominator
    
    # Take the ceiling to find the required number of iterations
    num_iterations = math.ceil(result)
    
    # Output the equation and the final result step-by-step
    print("The formula for the number of iterations (k) is:")
    print("k = ceil((c * T) / (2 * M))\n")
    
    print("Substituting the given values:")
    # Using f-strings to print the equation with the actual numbers
    print(f"k = ceil(({c} * {T}) / (2 * {M}))")
    print(f"k = ceil({numerator} / {denominator})")
    print(f"k = ceil({result:.4f})")
    print(f"k = {num_iterations}\n")

    print("-" * 70)
    print(f"The required number of iterations is: {num_iterations}")
    print(f"(The iteration counter will take values from 0 to {num_iterations - 1}).")
