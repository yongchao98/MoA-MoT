import numpy as np

# --- 1. Define Network Parameters ---
# Number of neurons in each input population (N_v, N_s)
N = 100

# Average firing probability (mu) for an input neuron.
# Derived from a Poisson process with an inter-event interval of 150s,
# and a discrete time step of 1ms (0.001s).
# The rate lambda = 1 event / 150 s.
# mu = probability of event in one time step = lambda * dt
dt = 0.001  # seconds
inter_event_interval = 150  # seconds
mu = dt / inter_event_interval

# Plasticity threshold (theta)
theta = 0.01

print(f"Calculating the critical covariance C based on the following parameters:")
print(f"Number of neurons N = {N}")
print(f"Average firing probability mu = {mu:.10f}")
print(f"Plasticity threshold theta = {theta}")
print("-" * 30)

# --- 2. Calculate Critical Covariance ---
# The formula derived from mean-field analysis is:
# C = mu * (2*N*theta - 1 - (2*N - 1)*mu)
# This C represents the covariance Cov(v_i, s_i) needed to balance
# potentiation and depression across the network.

term1 = 2 * N * theta
term2 = 1
term3 = (2 * N - 1) * mu

critical_covariance = mu * (term1 - term2 - term3)

# --- 3. Output the Final Equation and Result ---
print("The balancing condition leads to the equation for the critical covariance C:")
print("C = μ * (2*N*θ - 1 - (2*N - 1)*μ)\n")

print("Substituting the numerical values into the equation:")
# Using f-string formatting to display the full equation with numbers
final_equation_str = (
    f"C = {mu:.10f} * (2 * {N} * {theta} - 1 - (2 * {N} - 1) * {mu:.10f})"
)
print(final_equation_str)

# Show the intermediate calculation
print(f"C = {mu:.10f} * ({term1} - {term2} - {term3:.10f})")
print(f"C = {mu:.10f} * ({term1 - term2 - term3:.10f})")

print("\nFinal Result:")
print(f"The critical amount of correlation (covariance) is: C = {critical_covariance:.10f}")

# The problem asks for the final answer in a specific format
final_answer = critical_covariance
# <<<f"{final_answer}">>>