import numpy as np

# Set a seed for reproducibility
np.random.seed(42)

# --- Parameters ---
# Sample size at each step (must be > 1 for unbiased variance)
n = 10
# Number of steps to simulate and display
T = 5

# --- Theoretical Analysis ---
print("Analyzing the convergence of S_t = sum_{i=0 to t} sigma_i^2")
print("----------------------------------------------------------")
print("1. Convergence in L1 (Mean Absolute Error):")
print("Let sigma_t^2 be the variance at step t. By definition, it's an unbiased estimator of sigma_{t-1}^2.")
print("This implies E[sigma_t^2 | sigma_{t-1}^2] = sigma_{t-1}^2.")
print("Taking the full expectation, E[sigma_t^2] = E[E[sigma_t^2 | sigma_{t-1}^2]] = E[sigma_{t-1}^2].")
print(f"By induction, E[sigma_t^2] = E[sigma_0^2] = 1 for all t.")
print("The expectation of S_t is E[S_t] = E[sum_{i=0 to t} sigma_i^2] = sum_{i=0 to t} E[sigma_i^2] = t + 1.")
print("As t -> infinity, E[S_t] -> infinity. If S_t converged in L1 to a limit S, E[S_t] would need to converge to E[S].")
print("Since E[S_t] diverges, S_t does not converge in L1.")

print("\n2. Convergence in Distribution:")
print("S_t is a sum of non-negative terms (sigma_i^2 >= 0). This means the sequence S_t is non-decreasing.")
print("A non-decreasing sequence of random variables must converge almost surely to a limit S_inf (which could be +infinity).")
print("By the Monotone Convergence Theorem, E[S_inf] = E[lim S_t] = lim E[S_t] = lim (t + 1) = infinity.")
print("A random variable with an infinite expectation must be infinite with non-zero probability. For a non-decreasing sequence, this implies S_t diverges to infinity almost surely.")
print("Therefore, S_t does not converge in distribution to a proper (i.e., finite) random variable.")
print("----------------------------------------------------------\n")

print(f"Conclusion: S_t converges neither in L1 nor in distribution. Let's demonstrate with a simulation for T = {T} steps and n = {n}.\n")

# --- Simulation ---
mu = 0.0
sigma2 = 1.0

# Store the history of variances
sigma2_history = [sigma2]

# Run simulation
for t in range(1, T + 1):
    # Ensure variance is positive before sampling
    current_std_dev = np.sqrt(max(0, sigma2))
    
    # 1. Sample n variables from N(mu, sigma2)
    samples = np.random.normal(loc=mu, scale=current_std_dev, size=n)
    
    # 2. Compute new mu and unbiased sigma2
    mu = np.mean(samples)
    # Use ddof=1 for the unbiased estimator S^2 = 1/(n-1) * sum(x_i - x_bar)^2
    sigma2 = np.var(samples, ddof=1)
    
    sigma2_history.append(sigma2)

# --- Output the "final equation" ---
print(f"The simulated values of sigma_i^2 for this run are:")
for i, s2 in enumerate(sigma2_history):
    print(f"σ_{i}^2 = {s2:.4f}")

print(f"\nThe sum S_t at t = {T} is constructed as follows:")
equation_parts = [f"{s2:.4f}" for s2 in sigma2_history]
total_sum = sum(sigma2_history)
print(f"S_{T} = {' + '.join(equation_parts)} = {total_sum:.4f}")