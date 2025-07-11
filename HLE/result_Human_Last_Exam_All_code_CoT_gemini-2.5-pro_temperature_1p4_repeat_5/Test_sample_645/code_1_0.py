import numpy as np
import sympy
import math

def analyze_sets():
    """
    Generates terms and ratios for the four sets to illustrate their properties.
    """
    num_terms = 20

    # --- Set 1: Sums of Poisson random variables ---
    np.random.seed(0) # for reproducibility
    poisson_vars = np.random.poisson(1, size=num_terms * 2) # Generate more to ensure enough unique sums
    s1_list = sorted(list(set(np.cumsum(poisson_vars))))[:num_terms]
    ratios1 = [s1_list[i+1] / s1_list[i] for i in range(len(s1_list)-1)]

    print("--- Set 1: {sum(N_k)} where N_k ~ Poi(1) ---")
    print("Set (first {} terms):".format(num_terms))
    print(s1_list)
    print("Ratios s_{n+1}/s_n:")
    print([round(r, 2) for r in ratios1])
    print("The ratio appears to approach 1.\n")

    # --- Set 2: Powers of integers ---
    k = 4
    s2_list = [n**k for n in range(1, num_terms + 1)]
    ratios2 = [s2_list[i+1] / s2_list[i] for i in range(len(s2_list)-1)]
    
    print("--- Set 2: {n^k} for k=4 ---")
    print("Set (first {} terms):".format(num_terms))
    print(s2_list)
    print("Ratios s_{n+1}/s_n:")
    print([round(r, 2) for r in ratios2])
    print("The ratio appears to approach 1.\n")

    # --- Set 3: Prime numbers ---
    s3_list = list(sympy.primerange(1, sympy.prime(num_terms + 1)))
    ratios3 = [s3_list[i+1] / s3_list[i] for i in range(len(s3_list)-1)]

    print("--- Set 3: The set of primes ---")
    print("Set (first {} terms):".format(num_terms))
    print(s3_list)
    print("Ratios s_{n+1}/s_n:")
    print([round(r, 2) for r in ratios3])
    print("The ratio appears to approach 1.\n")

    # --- Set 4: floor((pi/2)^n) ---
    q = math.pi / 2
    s4_list = []
    # Using a while loop to ensure we get num_terms unique values
    n = 1
    while len(s4_list) < num_terms:
        val = math.floor(q**n)
        if not s4_list or val > s4_list[-1]:
             s4_list.append(val)
        n += 1

    ratios4 = [s4_list[i+1] / s4_list[i] for i in range(len(s4_list)-1) if s4_list[i] != 0]

    print("--- Set 4: {floor((pi/2)^n)} ---")
    print("Set (first {} terms):".format(num_terms))
    print(s4_list)
    print("Ratios s_{n+1}/s_n:")
    print([round(r, 2) for r in ratios4])
    print("The ratio appears to approach pi/2 = {:.2f} > 1.\n".format(q))
    
    print("Based on the analysis, sets 1, 2, and 3 have the required property, while set 4 does not.")

analyze_sets()