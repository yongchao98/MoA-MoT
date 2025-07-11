import numpy as np

def simulate_process_and_get_final_sum(n, T):
    """
    Simulates the iterative process for T steps and returns the final sum S_T.
    
    Args:
        n (int): Sample size at each step.
        T (int): Number of steps.
        
    Returns:
        float: The value of S_T = sum_{i=0 to T} sigma_i^2.
    """
    mu_current = 0.0
    sigma2_current = 1.0
    
    # Start the sum with the initial sigma_0^2
    S_t = sigma2_current
    
    for t in range(1, T + 1):
        # np.random.normal takes std dev, so we need sqrt(sigma^2).
        # Add a small epsilon for numerical stability if sigma2_current becomes zero.
        scale = np.sqrt(max(sigma2_current, 1e-9))
        samples = np.random.normal(loc=mu_current, scale=scale, size=n)
        
        # Compute new mu_t and sigma_t^2
        mu_next = np.mean(samples)
        # ddof=1 ensures an unbiased estimator for the variance
        sigma2_next = np.var(samples, ddof=1)
        
        # Update for the next iteration
        mu_current = mu_next
        sigma2_current = sigma2_next
        
        # Add the new variance to the sum
        S_t += sigma2_current
        
    return S_t

# --- Parameters ---
n = 10                  # Sample size, must be > 1
T = 5                   # Number of steps
num_simulations = 20000 # Number of simulations to run for averaging

# --- Simulation ---
final_S_T_values = []
for i in range(num_simulations):
    final_sum = simulate_process_and_get_final_sum(n, T)
    final_S_T_values.append(final_sum)

# --- Analysis and Output ---
simulated_mean_S_T = np.mean(final_S_T_values)
theoretical_mean_S_T = float(T + 1)

print(f"Simulation Analysis for n={n} and T={T}")
print("-" * 50)

# Print the final equation for the theoretical expectation
print("Theoretical Expectation Calculation:")
# E[S_T] = E[σ_0²] + E[σ_1²] + ... + E[σ_T²]
# Since σ_0²=1 and E[σ_i²]=1 for i>0, the sum is 1+1+...+1 (T+1 times)
expectation_terms = [1.0] * (T + 1)
equation_str = " + ".join(map(str, expectation_terms))
print(f"E[S_{T}] = E[σ_0²] + E[σ_1²] + ... + E[σ_{T}²]")
print(f"E[S_{T}] = {equation_str}")
print(f"E[S_{T}] = {theoretical_mean_S_T}")

print("-" * 50)
print(f"Comparison:")
print(f"Theoretical E[S_T] = {theoretical_mean_S_T}")
print(f"Simulated Average S_T after {num_simulations} runs = {simulated_mean_S_T:.4f}")

# Let's define the convergence options:
# A: converges in distribution only
# B: converges in L1 only
# C: converges in both
# D: converges in neither
