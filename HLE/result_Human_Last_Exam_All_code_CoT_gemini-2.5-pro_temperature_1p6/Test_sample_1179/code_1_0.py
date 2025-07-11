import numpy as np

# Set a seed for reproducibility
np.random.seed(42)

# Simulation parameters
n = 10          # Sample size at each step
T_max = 100     # Maximum number of time steps
num_paths = 10000 # Number of simulation paths to run

# --- Initialization ---
# mu_t for each path. mu_0 = 0
mu = np.zeros(num_paths)
# sigma2_t for each path. sigma_0^2 = 1
sigma2 = np.ones(num_paths)
# S_t = sum of sigma_i^2 from i=0 to t. S_0 = sigma_0^2 = 1
S = np.ones(num_paths)

# Store results for printing at these specific time steps
print_steps = [10, 20, 50, 100]

print(f"Running simulation with n={n} for {T_max} steps over {num_paths} paths.")
print("The analysis shows E[S_t] = t + 1. We compare this with the simulation average.")
print("-" * 75)
theoretical_S0 = 0 + 1
print(f"t=0:   E[S_t] = 0 + 1 = {theoretical_S0:<5.2f} (theoretical), {np.mean(S):.4f} (simulated)")

# --- Simulation Loop ---
for t in range(1, T_max + 1):
    # For each path, generate n samples from N(mu, sigma2)
    # np.random.normal takes standard deviation (scale), so we use sqrt(sigma2)
    samples = np.random.normal(loc=mu[:, np.newaxis], scale=np.sqrt(sigma2)[:, np.newaxis], size=(num_paths, n))

    # Update mu_t (sample mean) and sigma2_t (unbiased sample variance)
    mu = np.mean(samples, axis=1)
    sigma2 = np.var(samples, axis=1, ddof=1)
    
    # In rare cases, if all n samples are identical, variance is 0.
    # It cannot be negative.
    sigma2[sigma2 < 0] = 0

    # Update the sum S_t by adding the new variance
    S += sigma2

    # --- Output Results ---
    # Print results at the specified steps
    if t in print_steps:
        theoretical_E_S = t + 1
        simulated_E_S = np.mean(S)
        print(f"t={t:<3}: E[S_t] = {t} + 1 = {theoretical_E_S:<5.2f} (theoretical), {simulated_E_S:.4f} (simulated)")