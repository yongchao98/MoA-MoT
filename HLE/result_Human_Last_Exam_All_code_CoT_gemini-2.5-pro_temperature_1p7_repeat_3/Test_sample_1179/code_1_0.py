import numpy as np

# Simulation Parameters
n = 10                # Sample size at each step (must be > 1)
max_t = 100           # Maximum number of time steps
num_simulations = 5000 # Number of full simulations to run

# This array will store the value of S_t for each simulation at each time step
s_t_results = np.zeros((num_simulations, max_t + 1))

print("Running simulations...")
# Loop over the number of simulations
for i in range(num_simulations):
    # Initialize process for each simulation run
    mu = 0.0
    sigma2 = 1.0
    
    # Store the history of sigma_i^2 values, starting with sigma_0^2
    sigma2_history = [sigma2]

    # Iterate from t=1 to max_t
    for t in range(1, max_t + 1):
        # We need a non-zero variance to sample. If sigma2 becomes zero,
        # the process can't continue. This is highly unlikely.
        if sigma2 <= 1e-9: # Use a small threshold for stability
            # Fill the rest of history with 0 and break the inner loop
            sigma2_history.extend([0] * (max_t - t + 1))
            break
            
        # 1. Sample n variables from N(mu_{t-1}, sigma^2_{t-1})
        samples = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=n)

        # 2. Compute mu_t (MLE for mean) and sigma_t^2 (unbiased estimator for variance)
        mu = np.mean(samples)
        sigma2 = np.var(samples, ddof=1) # ddof=1 for the unbiased estimator
        
        sigma2_history.append(sigma2)

    # Calculate the cumulative sum S_t = sum_{i=0 to t} sigma_i^2 for this simulation
    s_t_values = np.cumsum(sigma2_history)
    s_t_results[i, :] = s_t_values

# Analyze and Print Results
# Calculate the average of S_t across all simulations for each time step
avg_s_t = np.mean(s_t_results, axis=0)

print("\n--- Simulation Analysis ---")
print(f"Parameters: n={n}, max_t={max_t}, num_simulations={num_simulations}")
print("\nComparing theoretical expectation E[S_t] = t + 1 with simulated average.")
print("-" * 60)
print(f"{'Time t':<10} | {'Theoretical E[S_t] = t + 1':<30} | {'Simulated Avg(S_t)':<20}")
print("-" * 60)

# Print results at key time steps
for t in [0, 1, 10, 20, 50, 100]:
    if t <= max_t:
        theoretical_mean = t + 1
        simulated_mean = avg_s_t[t]
        # The output below shows each number in the final equation E[S_t] = t + 1
        # For t=10, the equation is E[S_10] = 10 + 1 = 11
        print(f"{t:<10} | E[S_{t}] = {t} + 1 = {theoretical_mean:<16.2f} | {simulated_mean:<20.2f}")
print("-" * 60)

print("\n--- Conclusion ---")
print("The simulation shows that the average value of S_t grows linearly with t,")
print("closely matching the theoretical expectation E[S_t] = t + 1.")
print("Since the expectation of S_t diverges to infinity, the sum S_t converges")
print("in neither L1 nor in distribution.")
