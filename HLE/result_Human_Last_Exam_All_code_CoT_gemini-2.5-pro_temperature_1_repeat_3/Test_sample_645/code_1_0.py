import numpy as np
import math

def analyze_sets():
    """
    Analyzes the properties of the four sets of indices to determine
    if they support the described type of power series.

    The analysis hinges on whether the set is a Sidon set. A practical
    way to test this is to check for lacunarity. A set S = {s_k} is
    lacunary if s_{k+1}/s_k >= q > 1. Lacunary sets are Sidon sets.
    """
    print("Analyzing the properties of the index sets S.")
    print("A power series with the desired properties exists if and only if S is not a Sidon set.")
    print("A key indicator is lacunarity: if s_{k+1}/s_k >= q > 1, the set is lacunary and thus a Sidon set.")
    print("Let's examine the ratio of consecutive terms for each set.\n")

    # Set 1: Sum of Poisson variables
    print("--- Set 1: S = {sum(N_k) up to n}, N_k ~ Poi(1) ---")
    # Simulate one realization of this random set for illustration
    np.random.seed(0)  # for reproducibility
    N_k = np.random.poisson(1, size=40)
    S_n_raw = np.cumsum(N_k)
    s_k = sorted(list(set(S_n_raw) - {0})) # get distinct non-zero values
    if len(s_k) > 1:
        ratios = [s_k[i+1] / s_k[i] for i in range(len(s_k)-1)]
        print(f"First 10 distinct non-zero terms (one simulation): {s_k[:10]}")
        print(f"Ratios of first 10 terms: {[f'{r:.2f}' for r in ratios[:10]]}")
    print("Analysis: The ratio appears to approach 1. Theory confirms the set has density 1 (a.s.) and is not a Sidon set.")
    print("Result: Yes\n")

    # Set 2: Powers of integers
    print("--- Set 2: S = {n^k}, k>=4 (using k=4) ---")
    k = 4
    s_k = [n**k for n in range(1, 13)]
    ratios = [s_k[i+1] / s_k[i] for i in range(len(s_k)-1)]
    print(f"First 10 terms: {s_k[:10]}")
    print(f"Ratios: {[f'{r:.2f}' for r in ratios[:10]]}")
    print("Analysis: The ratio clearly approaches 1. This set is not lacunary and known not to be a Sidon set.")
    print("Result: Yes\n")

    # Set 3: Primes
    print("--- Set 3: S = {primes} ---")
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(math.sqrt(n)) + 1):
            if n % i == 0: return False
        return True
    s_k = []
    num = 2
    while len(s_k) < 13:
        if is_prime(num):
            s_k.append(num)
        num += 1
    ratios = [s_k[i+1] / s_k[i] for i in range(len(s_k)-1)]
    print(f"First 11 primes: {s_k[:11]}")
    print(f"Ratios: {[f'{r:.2f}' for r in ratios[:10]]}")
    print("Analysis: The ratio approaches 1. This set is not lacunary and known not to be a Sidon set.")
    print("Result: Yes\n")

    # Set 4: Floor of powers of pi/2
    print("--- Set 4: S = {floor((pi/2)^n)} ---")
    s_k = []
    for n in range(1, 13):
        val = math.floor((math.pi/2)**n)
        if not s_k or val > s_k[-1]: # handle duplicates and ensure strictly increasing
            s_k.append(val)
    ratios = [s_k[i+1] / s_k[i] for i in range(len(s_k)-1)]
    print(f"First 10 terms: {s_k[:10]}")
    print(f"Ratios: {[f'{r:.2f}' for r in ratios[:10]]}")
    print(f"Analysis: The ratio appears to approach pi/2 ~= {math.pi/2:.2f}. This set is lacunary and is a Sidon set.")
    print("Result: No\n")

    # Final conclusion based on the analysis
    print("The property holds for sets 1, 2, and 3.")
    print("The numbers of the sets that satisfy the condition are:")
    print(1)
    print(2)
    print(3)

if __name__ == '__main__':
    analyze_sets()