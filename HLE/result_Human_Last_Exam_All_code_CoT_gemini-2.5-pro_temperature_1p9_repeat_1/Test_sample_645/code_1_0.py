import numpy as np
import math

def analyze_sets():
    """
    Analyzes four sets of integers to determine if they can be the index set 
    for a power series that converges on the closed unit disk but not absolutely.

    The core mathematical property is that such a power series exists if and only if
    the set of indices S is NOT a "Sidon set". A practical way to test this is
    to check if the set is "lacunary" (has large gaps). A set S={s_n} is 
    lacunary if s_{n+1}/s_n >= q > 1. Lacunary sets are always Sidon sets.

    This script checks each set for the lacunary property to support the conclusion.
    """

    def check_lacunarity(s, name):
        """
        Analyzes a set S to see if it is lacunary by computing the ratio s_{n+1}/s_n.
        """
        print(f"--- Analyzing Set: {name} ---")
        # Ensure the set is sorted, contains unique elements, and is positive
        s_unique_sorted = sorted([val for val in list(set(s)) if val > 0])
        
        if len(s_unique_sorted) < 2:
            print("Set has fewer than 2 elements, cannot check ratios.")
            return False

        ratios = [s_unique_sorted[i+1] / s_unique_sorted[i] for i in range(len(s_unique_sorted) - 1)]
        
        print(f"First 10 terms of the set: {s_unique_sorted[:10]}")
        print(f"Ratios s_(n+1)/s_n for first 10 terms: {[f'{r:.3f}' for r in ratios[:10]]}")
        
        # Check the behavior of the ratio for the last few computed terms
        if len(ratios) > 20:
             print(f"Last 5 computed ratios: {[f'{r:.3f}' for r in ratios[-5:]]}")
        
        last_ratio = ratios[-1]
        if last_ratio > 1.1:
            print(f"Conclusion: The ratio appears to approach a limit > 1 ({last_ratio:.3f}). The set is LACUNARY.")
            print("This implies it IS a Sidon set. The desired power series does NOT exist.\n")
            return True
        else:
            print(f"Conclusion: The ratio appears to approach 1 ({last_ratio:.3f}). The set is NOT lacunary.")
            print("This suggests it is NOT a Sidon set. The desired power series EXISTS.\n")
            return False

    # 1. S = { sum_{k<=n} N_k } where N_k ~ Poi(1)
    np.random.seed(42) # for reproducibility
    N_terms = 500
    N_k = np.random.poisson(lam=1.0, size=N_terms)
    S1_raw = np.cumsum(N_k)
    check_lacunarity(S1_raw, "1. S = {sums of Poisson(1) variables}")

    # 2. S = { n^k } for k >= 4
    k = 4
    n_values = np.arange(1, 101)
    S2 = [n**k for n in n_values]
    check_lacunarity(S2, f"2. S = {{n^k}} for k={k}")

    # 3. The set of primes
    def get_primes(n_max):
        sieve = np.ones(n_max + 1, dtype=bool)
        sieve[0:2] = False
        for p in range(2, int(np.sqrt(n_max)) + 1):
            if sieve[p]:
                sieve[p*p::p] = False
        return list(np.where(sieve)[0])
    
    S3 = get_primes(1000)
    check_lacunarity(S3, "3. The set of primes")
    
    # 4. S = { floor( (pi/2)^n ) }
    N_terms_4 = 40
    S4 = [math.floor((math.pi/2)**n) for n in range(1, N_terms_4 + 1)]
    check_lacunarity(S4, "4. S = {floor((pi/2)^n)}")

    print("Summary:")
    print("Sets 1, 2, and 3 are NOT lacunary and are known not to be Sidon sets.")
    print("Set 4 IS lacunary and therefore is a Sidon set.")
    print("A power series with the desired properties exists if and only if the index set is NOT a Sidon set.")
    print("Therefore, the property holds for sets 1, 2, and 3 only.")

# Execute the analysis
analyze_sets()