import math
import numpy as np

def analyze_sets():
    """
    Generates terms and computes ratios for four different sets of integers
    to analyze their properties, particularly lacunarity.
    """
    n_terms = 30 # Number of terms to generate for each set

    # --- Set 1: Sums of Poisson random variables ---
    print("--- Set 1: S = {sum(N_k) for N_k ~ Poi(1)} ---")
    np.random.seed(0) # for reproducibility
    current_sum = 0
    set1 = []
    attempts = 0
    # Generate unique, positive terms
    while len(set1) < n_terms and attempts < 1000:
        current_sum += np.random.poisson(1)
        if current_sum > 0 and (len(set1) == 0 or current_sum > set1[-1]):
            set1.append(current_sum)
        attempts += 1
    
    print("First {} terms:".format(len(set1)))
    print(set1)
    
    # Calculate ratios
    ratios1 = [set1[i] / set1[i-1] for i in range(1, len(set1))]
    print("Ratios of consecutive terms s_n / s_{n-1}:")
    # Print ratios in a formatted way
    print(" ".join(["{:.2f}".format(r) for r in ratios1]))
    print("Observation: The ratios tend towards 1, so the set is not lacunary.\n")

    # --- Set 2: k-th powers (k=4) ---
    print("--- Set 2: S = {n^4} ---")
    set2 = [n**4 for n in range(1, n_terms + 1)]
    print("First {} terms:".format(len(set2)))
    print(set2)
    
    # Calculate ratios
    ratios2 = [set2[i] / set2[i-1] for i in range(1, len(set2))]
    print("Ratios of consecutive terms s_n / s_{n-1}:")
    print(" ".join(["{:.2f}".format(r) for r in ratios2]))
    print("Observation: The ratios tend towards 1, so the set is not lacunary.\n")

    # --- Set 3: Prime numbers ---
    print("--- Set 3: The set of primes ---")
    def is_prime(num):
        if num < 2:
            return False
        for i in range(2, int(math.sqrt(num)) + 1):
            if num % i == 0:
                return False
        return True

    primes = []
    num = 2
    while len(primes) < n_terms:
        if is_prime(num):
            primes.append(num)
        num += 1
    
    print("First {} terms:".format(len(primes)))
    print(primes)
    
    # Calculate ratios
    ratios3 = [primes[i] / primes[i-1] for i in range(1, len(primes))]
    print("Ratios of consecutive terms s_n / s_{n-1}:")
    print(" ".join(["{:.2f}".format(r) for r in ratios3]))
    print("Observation: The ratios tend towards 1, so the set is not lacunary.\n")

    # --- Set 4: Floor of powers of pi/2 ---
    print("--- Set 4: S = {floor((pi/2)^n)} ---")
    set4 = []
    val = 1
    # Use range that will likely produce enough unique terms
    for n in range(1, n_terms * 2):
        val = math.floor((math.pi/2)**n)
        if val > 0 and (len(set4) == 0 or val > set4[-1]):
             set4.append(val)
        if len(set4) == n_terms:
            break

    print("First {} terms:".format(len(set4)))
    print(set4)

    # Calculate ratios
    ratios4 = [set4[i] / set4[i-1] for i in range(1, len(set4))]
    print("Ratios of consecutive terms s_n / s_{n-1}:")
    print(" ".join(["{:.2f}".format(r) for r in ratios4]))
    print("Observation: The ratios approach pi/2 (~1.57), which is > 1. The set is lacunary.\n")

if __name__ == '__main__':
    analyze_sets()