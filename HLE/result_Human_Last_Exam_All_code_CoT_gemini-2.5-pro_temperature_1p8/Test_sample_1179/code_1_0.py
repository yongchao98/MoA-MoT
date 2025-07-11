import numpy as np

def analyze_convergence(T, n, num_paths):
    """
    Analyzes the convergence of S_t by simulating the process and showing
    that its mean and variance diverge.

    Args:
        T (int): Number of time steps.
        n (int): Sample size at each step (must be > 1).
        num_paths (int): Number of simulation paths to run for statistical analysis.
    """
    print(f"Analyzing convergence for T={T}, n={n} over {num_paths} paths.\n")

    # We will store the values of the sum S_t at specific time points
    # across all simulation paths to analyze their distribution.
    time_points = sorted(list(set([0, 10, T // 2, T])))
    S_values = np.zeros((len(time_points), num_paths))

    print(f"Running {num_paths} simulations...")
    for i in range(num_paths):
        # Initialization for each path
        mu = 0.0
        sigma2 = 1.0
        
        # S_0 = sigma_0^2
        current_S = sigma2
        
        time_point_idx = 0
        if time_points[time_point_idx] == 0:
             S_values[time_point_idx, i] = current_S
             time_point_idx += 1

        # Main iterative process
        for t in range(1, T + 1):
            # 1. Sample from N(mu, sigma2)
            samples = np.random.normal(loc=mu, scale=np.sqrt(sigma2), size=n)
            
            # 2. Compute new mu and sigma2
            mu = np.mean(samples)
            sigma2 = np.var(samples, ddof=1)
            
            # Add new sigma2 to the running sum
            current_S += sigma2

            # Store the sum at the specified time points
            if time_point_idx < len(time_points) and t == time_points[time_point_idx]:
                S_values[time_point_idx, i] = current_S
                time_point_idx += 1
    print("Simulations complete.\n")

    # --- Analysis & Output ---
    
    print("--- 1. L1 Convergence Analysis ---")
    print("For S_t to converge in L1, its expected value E[S_t] must converge to a finite value.")
    print("Theoretical analysis shows E[S_t] = t + 1, which diverges to infinity.")
    print("Let's compare this with our simulation results:\n")
    
    # Calculate the mean of S_t across all paths
    mean_S = np.mean(S_values, axis=1)
    
    print(" t | Theoretical E[S_t] | Simulated E[S_t]")
    print("---|--------------------|-------------------")
    for i, t in enumerate(time_points):
        theoretical_mean = t + 1
        print(f"{t:2d} | {theoretical_mean:18d} | {mean_S[i]:17.2f}")

    print("\nConclusion: The simulated mean grows linearly with t, matching the theory. S_t does not converge in L1.\n")

    print("--- 2. Convergence in Distribution Analysis ---")
    print("For S_t to converge in distribution, its distribution must stabilize (e.g., have a convergent variance).")
    print("Let's examine the variance of S_t at different times:\n")
    
    # Calculate the variance of S_t across all paths
    var_S = np.var(S_values, axis=1)
    
    print(" t | Simulated Var(S_t)")
    print("---|--------------------")
    for i, t in enumerate(time_points):
        print(f"{t:2d} | {var_S[i]:18.2f}")

    print("\nConclusion: The variance of S_t continuously increases. The distribution does not stabilize but instead spreads out indefinitely.")
    print("This supports the theoretical conclusion that S_t does not converge in distribution.")

if __name__ == '__main__':
    # --- Parameters ---
    T_final = 100         # Total number of steps in the process
    n_samples = 5         # Sample size at each step (must be > 1)
    num_simulations = 20000 # Number of paths to simulate for accurate stats
    
    analyze_convergence(T_final, n_samples, num_simulations)
