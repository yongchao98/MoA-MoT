import numpy as np
import math

def is_prime(num):
    """Checks if a number is prime."""
    if num < 2:
        return False
    for i in range(2, int(math.sqrt(num)) + 1):
        if num % i == 0:
            return False
    return True

def generate_set_1(num_terms):
    """
    S = { sum_{k <= n} N_k : n in N } where N_k ~ Poi(1).
    We generate one realization of this random process.
    """
    np.random.seed(42) # for reproducibility
    S = []
    current_sum = 0
    # Use while loop to ensure we get num_terms distinct values
    n = 0
    while len(S) < num_terms:
        # To avoid zero sums, start with a non-zero value
        if n == 0:
            current_sum = np.random.poisson(1)
            if current_sum == 0: continue
        else:
            current_sum += np.random.poisson(1)
        
        if current_sum not in S: # make sure elements are unique
             S.append(current_sum)
        n += 1
    return S

def generate_set_2(num_terms, k=4):
    """
    S = { (n^k) : n in N} for k >= 4.
    """
    return [n**k for n in range(1, num_terms + 1)]

def generate_set_3(num_terms):
    """
    S = the set of primes.
    """
    S = []
    num = 2
    while len(S) < num_terms:
        if is_prime(num):
            S.append(num)
        num += 1
    return S

def generate_set_4(num_terms):
    """
    S = { floor( (pi/2)^n ) : n in N }.
    """
    S = []
    n = 1
    # Use while loop to ensure we get num_terms distinct values
    while len(S) < num_terms:
      val = math.floor((math.pi / 2)**n)
      if val > 0 and (len(S) == 0 or val > S[-1]):
        S.append(val)
      n += 1
    return S

def analyze_set(name, s_list):
    """Analyzes and prints properties of a given set."""
    print(f"--- Analysis for Set {name} ---")
    if len(s_list) < 2:
        print("Set has fewer than 2 elements, cannot compute ratios.")
        print("\n")
        return
    
    print(f"First {min(len(s_list), 15)} elements: {s_list[:15]}")
    ratios = [s_list[i+1] / s_list[i] for i in range(len(s_list) - 1)]
    print(f"Ratios of consecutive terms (first {min(len(ratios), 10)}):")
    print([round(r, 2) for r in ratios[:10]])

    # Check the trend in the latter half of the ratios
    if len(ratios) > 10:
        later_ratios = ratios[len(ratios)//2:]
        if all(r > 1.2 for r in later_ratios): # Check if bounded away from 1
             print("The ratio s_{n+1}/s_n appears bounded below by a constant q > 1.")
             print("This suggests the set IS a Hadamard lacunary sequence and thus a Sidon set.")
        else:
             print("The ratio s_{n+1}/s_n appears to approach 1.")
             print("This suggests the set IS NOT a Hadamard lacunary sequence.")
    else:
        # For small number of terms, just check the last ratio
        if ratios[-1] > 1.2:
            print("The ratio s_{n+1}/s_n appears bounded below by a constant q > 1.")
            print("This suggests the set IS a Hadamard lacunary sequence and thus a Sidon set.")
        else:
            print("The ratio s_{n+1}/s_n appears to approach 1.")
            print("This suggests the set IS NOT a Hadamard lacunary sequence.")
    print("\n")


def solve():
    """Main function to perform analysis and print conclusions."""
    num_terms_std = 30
    num_terms_fast_growth = 20
    
    # Header
    print("A power series with the desired properties exists if and only if the set of indices S is NOT a Sidon set.")
    print("Hadamard lacunary sequences (where s_{n+1}/s_n >= q > 1) are Sidon sets.")
    print("Let's numerically investigate this lacunarity property for each set.\n")

    # Set 1
    s1 = generate_set_1(num_terms_std)
    analyze_set("1: Sums of Poisson random variables", s1)
    
    # Set 2
    s2 = generate_set_2(num_terms_std, k=4)
    analyze_set("2: Fourth powers of integers", s2)
    
    # Set 3
    s3 = generate_set_3(num_terms_std)
    analyze_set("3: Prime numbers", s3)
    
    # Set 4
    s4 = generate_set_4(num_terms_fast_growth)
    analyze_set("4: Floor of powers of pi/2", s4)
    
    # Final Conclusion
    print("--------------------------------------------------")
    print("CONCLUSION:")
    print("Sets 1, 2, and 3 are NOT Sidon sets based on known theorems.")
    print("Set 4 IS a Sidon set (since it is Hadamard lacunary).")
    print("\nTherefore, a power series with the specified convergence properties exists for sets 1, 2, and 3.")
    print("The correct option includes the numbers: 1, 2, 3")

solve()