import numpy as np

def analyze_variance_sum(n, num_steps):
    """
    Simulates the iterative process and analyzes the convergence of the sum of variances.

    Args:
        n (int): The number of samples to draw at each step. Must be > 1.
        num_steps (int): The total number of steps to simulate.
    """
    if n <= 1:
        print("Error: Sample size n must be greater than 1 for unbiased variance calculation.")
        return

    # Initial conditions
    mu_current = 0.0
    sigma2_current = 1.0

    # History tracking
    sigma2_history = [sigma2_current]
    
    # Iterative process
    for _ in range(num_steps):
        # Sample from N(mu_{t-1}, sigma^2_{t-1})
        # Note: np.random.normal takes standard deviation (sqrt of variance) as scale
        samples = np.random.normal(loc=mu_current, scale=np.sqrt(sigma2_current), size=n)

        # Compute MLE for mean and unbiased estimator for variance
        mu_current = np.mean(samples)
        # ddof=1 provides the unbiased estimator for variance
        sigma2_current = np.var(samples, ddof=1)
        
        # Guard against numerically unstable results (though unlikely here)
        if sigma2_current < 0:
            sigma2_current = 0
            
        sigma2_history.append(sigma2_current)

    # --- Output Results ---
    
    print(f"This script analyzes the convergence of S_t = sum(sigma_i^2) for i=0 to t.")
    print(f"Simulation run with n={n} samples for {num_steps} steps.")
    print("-" * 50)
    
    # Display the sum for the first 10 steps to satisfy the "show the equation" requirement.
    display_steps = min(10, num_steps)
    s_display = sum(sigma2_history[:display_steps+1])
    
    # Building the equation string
    equation_str = f"S_{display_steps} = "
    for i in range(display_steps + 1):
        equation_str += f"{sigma2_history[i]:.4f}"
        if i < display_steps:
            equation_str += " + "
    equation_str += f" = {s_display:.4f}"
    
    print("Example sum for the first few steps:")
    print(equation_str)
    print("-" * 50)

    # Calculate and print the final sum
    s_final = sum(sigma2_history)
    print(f"The total sum after {num_steps} steps is S_{num_steps} = {s_final:.4f}")
    print("-" * 50)

    # Print the final conclusion
    print("Conclusion on Convergence:")
    print("The sum S_t does NOT converge in L1 or in distribution.")
    print("\nReasoning:")
    print("1. The expected value of each term is E[sigma_i^2] = 1 for all i.")
    print("2. The expectation of the sum is E[S_t] = E[sum_{i=0 to t} sigma_i^2] = sum(E[sigma_i^2]) = t + 1.")
    print("3. As t -> infinity, E[S_t] -> infinity.")
    print("4. Convergence in L1 requires the expectation to converge to a finite value, which is not the case.")
    print("5. Convergence in distribution is also ruled out because the sum S_t grows unboundedly (diverges to infinity).")


if __name__ == '__main__':
    # --- Parameters ---
    # Number of samples per step (must be > 1)
    sample_size_n = 5
    # Total number of iterations
    total_steps_t = 1000

    analyze_variance_sum(n=sample_size_n, num_steps=total_steps_t)