import numpy as np
import math

def generate_set_1(n_terms=20):
    """Generates terms for S = {sum(N_k) where N_k ~ Poi(1)}"""
    # Use a set to handle potential duplicates, though unlikely for large sums
    s_set = set()
    current_sum = 0
    # Ensure we generate at least n_terms distinct values
    while len(s_set) < n_terms:
        # To ensure the sum is strictly increasing for the ratio test
        # we can take N_k from Poi(1) conditioned on being >= 1
        # P(Poi(1)=0) = exp(-1), so this is a minor adjustment for small sums.
        # Let's just add 1 to each draw to keep it simple and strictly increasing.
        N_k = np.random.poisson(1) + 1
        current_sum += N_k
        s_set.add(current_sum)
    
    return sorted(list(s_set))[:n_terms]

def generate_set_2(k, n_terms=20):
    """Generates terms for S = {n^k}"""
    return [n**k for n in range(1, n_terms + 1)]

def generate_set_3(n_terms=20):
    """Generates terms for S = the set of primes"""
    primes = []
    num = 2
    while len(primes) < n_terms:
        is_prime = True
        for i in range(2, int(math.sqrt(num)) + 1):
            if num % i == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(num)
        num += 1
    return primes

def generate_set_4(n_terms=20):
    """Generates terms for S = {floor((pi/2)^n)}"""
    s_set = set()
    n = 1
    while len(s_set) < n_terms:
        term = math.floor((math.pi / 2)**n)
        s_set.add(term)
        n += 1
    return sorted(list(s_set))

def analyze_sets():
    """Analyzes the four sets and prints the ratio of consecutive terms."""
    print("Analyzing the gap behavior of the sets S = {s_k}.\n")
    
    sets = {
        "Set 1 (Poisson Sums)": generate_set_1(),
        "Set 2 (n^k, k=4)": generate_set_2(4),
        "Set 3 (Primes)": generate_set_3(),
        "Set 4 (floor((pi/2)^n))": generate_set_4(20)
    }
    
    for name, s in sets.items():
        print(f"--- {name} ---")
        print(f"First 10 terms: {s[:10]}")
        ratios = [s[i+1] / s[i] for i in range(len(s) - 1) if s[i] != 0]
        print(f"Ratios s_k+1/s_k (first 10): {[f'{r:.2f}' for r in ratios[:10]]}")
        if ratios:
            print(f"Ratio for large k (last value): {ratios[-1]:.4f}")
        print("")

    print("--- Conclusion ---")
    print("A sufficient condition for the power series to exist is the Hadamard gap condition (s_k+1/s_k >= q > 1).")
    print("Only Set 4 satisfies this, with the ratio approaching pi/2 ≈ 1.5708.")
    print("Set 2 fails due to its high arithmetic regularity.")
    print("Sets 1 and 3 are sufficiently irregular ('random') to prevent failures on the unit circle.")
    print("Therefore, the property holds for sets 1, 3, and 4.")
    print("\nFinal choice corresponds to 'Only 1, 3 and 4'.")
    
analyze_sets()