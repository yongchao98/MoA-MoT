import math

# --- 1. Define constants based on the problem statement ---
# We are given epsilon and mu in terms of k_B*T, so we can work with the dimensionless
# quantities beta*epsilon and beta*mu, where beta = 1/(k_B*T).

# Given: epsilon = -k_B*T / (2*pi)
# Therefore, beta*epsilon = -1 / (2*pi)
beta_epsilon = -1.0 / (2.0 * math.pi)

# Given: mu = 0.1 * k_B*T
# Therefore, beta*mu = 0.1
beta_mu = 0.1

# The total coordination number z is the sum of horizontal and vertical neighbors.
z_horizontal = 4
z_vertical = 8
z = z_horizontal + z_vertical

# --- 2. Display the self-consistency equation with numerical values ---
print("The self-consistency equation for the average occupancy <n> is:")
print("<n> = 1 / [exp((epsilon * z * <n> - mu) / (k_B * T)) + 1]")
print("\nSubstituting beta*epsilon, beta*mu, and z with the given values:")
# We display the final equation with all the numbers that will be used.
# The number for beta_epsilon is calculated and shown truncated for readability.
print(f"<n> = 1 / [exp({beta_epsilon:.5f} * {z} * <n> - {beta_mu}) + 1]")

# --- 3. Solve the equation numerically using fixed-point iteration ---
# The equation to solve is: n = 1 / [exp(beta_epsilon * z * n - beta_mu) + 1]
# We define a function for the right-hand side.
def g(n_val):
    """Calculates the right-hand side of the self-consistency equation."""
    exponent = beta_epsilon * z * n_val - beta_mu
    return 1.0 / (math.exp(exponent) + 1.0)

# Initialize <n> with a starting guess.
# Since <n> represents an occupancy, it must be between 0 and 1. 0.5 is a common choice.
n_avg = 0.5

# Iterate until the value of n_avg converges to a stable solution.
# We iterate a fixed number of times, which is sufficient for convergence here.
for _ in range(100):
    n_avg = g(n_avg)

# --- 4. Print the final result ---
# The result is rounded to three decimal places as requested.
print("\nThe converged average occupancy per site <n> is:")
print(f"{n_avg:.3f}")

# The final answer in the required format
final_answer = round(n_avg, 3)
# print(f"\n<<<{final_answer}>>>") # This is for the final answer extraction