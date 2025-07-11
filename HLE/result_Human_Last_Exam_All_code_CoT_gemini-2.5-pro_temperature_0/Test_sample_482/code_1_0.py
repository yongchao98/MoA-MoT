import math

# --- Parameters ---
# These are example values since they were not specified in the problem.

# N_v: Number of neurons in input layer v
N_v = 100
# N_s: Number of neurons in input layer s
N_s = 100
# Total number of input neurons
N_total = N_v + N_s

# Inter-event interval for the Poisson process (in seconds)
inter_event_interval = 150.0
# The average rate 'mu' for a single neuron is the total rate (1/interval)
# divided by the total number of input neurons.
mu = (1.0 / inter-event_interval) / N_total

# theta: Heterosynaptic offset constant. A common choice for stability in
# similar models is 1 divided by the total number of inputs driving plasticity.
theta = 1.0 / N_total

# --- Calculation ---
# We use the derived formula for the critical covariance C_crit:
# C_crit = -mu * [1 + (N_v + N_s) * (mu - theta)]

# Calculate the term (mu - theta)
mu_minus_theta = mu - theta

# Calculate the term (N_v + N_s) * (mu - theta)
N_times_mu_minus_theta = N_total * mu_minus_theta

# Calculate the full expression
C_crit = -mu * (1 + N_times_mu_minus_theta)

# --- Output ---
print("### Determining the Critical Correlation ###\n")
print("This script calculates the critical covariance 'C' between input populations")
print("that balances potentiation and depression in the network.\n")

print("--- Model Parameters (Examples) ---")
print(f"Number of neurons in v (N_v): {N_v}")
print(f"Number of neurons in s (N_s): {N_s}")
print(f"Average input rate (mu): {mu:.6f}")
print(f"Heterosynaptic offset (theta): {theta:.6f}\n")

print("--- Calculation Steps ---")
print("The formula for the critical covariance is:")
print("C = -μ * [1 + (N_v + N_s) * (μ - θ)]\n")

print("Plugging in the values:")
# Using f-strings with specified precision for cleaner output
print(f"C = -{mu:.6f} * [1 + ({N_v} + {N_s}) * ({mu:.6f} - {theta:.6f})]")
print(f"C = -{mu:.6f} * [1 + {N_total} * ({mu_minus_theta:.6f})]")
print(f"C = -{mu:.6f} * [1 + {N_times_mu_minus_theta:.6f}]")
print(f"C = -{mu:.6f} * [{(1 + N_times_mu_minus_theta):.6f}]\n")

print("--- Final Result ---")
print(f"The critical covariance C is: {C_crit}")

# The correlation coefficient rho is C / (sigma_v * sigma_s).
# For Poisson-like statistics, sigma^2 = mu, so rho = C / mu.
rho_crit = C_crit / mu
print(f"The critical correlation coefficient rho is: {rho_crit:.6f}")
<<< -2.2222222222222223e-07>>>