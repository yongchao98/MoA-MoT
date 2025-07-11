import math
import numpy as np

def find_max_ap_length(s_list):
    """
    Finds the maximum length of an Arithmetic Progression in a sorted list of unique integers.
    This is an O(N^2) implementation, efficient enough for the demonstration.
    """
    if len(s_list) < 2:
        return len(s_list)
    
    s_set = set(s_list)
    max_len = 1
    if s_list:
        max_len = 2 if len(s_list) > 1 else 1

    for i in range(len(s_list)):
        for j in range(i + 1, len(s_list)):
            a1 = s_list[i]
            a2 = s_list[j]
            d = a2 - a1

            # Optimization: if a1-d exists, we already checked this AP
            if (a1 - d) in s_set:
                continue

            # Count the length of this AP
            current_len = 2
            next_term = a2 + d
            while next_term in s_set:
                current_len += 1
                next_term += d
            
            if current_len > max_len:
                max_len = current_len
    return max_len

def get_primes(n):
    """Generate primes up to n using a sieve."""
    primes = []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(math.sqrt(n)) + 1):
        if is_prime[p]:
            for multiple in range(p*p, n + 1, p):
                is_prime[multiple] = False
    for p in range(2, n + 1):
        if is_prime[p]:
            primes.append(p)
    return primes

print("Analyzing sets for the maximum length of arithmetic progressions (APs).\n")

# Case 1: Partial sums of Poisson variables
print("--- Case 1: Partial sums of Poisson variables ---")
np.random.seed(42) # for reproducibility
num_terms_1 = 5000
poisson_steps = np.random.poisson(1, num_terms_1)
s1_full = np.cumsum(poisson_steps)
s1 = sorted(list(set(s1_full)))
max_l_1 = find_max_ap_length(s1)
print(f"Generated {len(s1)} unique values up to {s1[-1]}.")
print(f"Found maximum AP length in sample: {max_l_1}")
print("CONCLUSION: Theory states this set contains arbitrarily long APs. The property does NOT hold.")
print("-" * 50)

# Case 2: k-th powers (k=4)
print("--- Case 2: k-th powers (k=4) ---")
k = 4
num_terms_2 = 100
s2 = [n**k for n in range(1, num_terms_2 + 1)]
max_l_2 = find_max_ap_length(s2)
print(f"Generated {len(s2)} values (1^{k} to {num_terms_2}^{k}).")
print(f"Found maximum AP length: {max_l_2}")
print("CONCLUSION: Theory states the max AP length is 2. The property HOLDS.")
print("-" * 50)

# Case 3: Prime numbers
print("--- Case 3: Prime numbers ---")
max_prime = 3000
s3 = get_primes(max_prime)
max_l_3 = find_max_ap_length(s3)
print(f"Generated {len(s3)} primes up to {max_prime}.")
print(f"Found maximum AP length in sample: {max_l_3}")
print("CONCLUSION: Theory states this set contains arbitrarily long APs. The property does NOT hold.")
print("-" * 50)

# Case 4: floor((pi/2)^n)
print("--- Case 4: floor((pi/2)^n) ---")
num_terms_4 = 40 
s4_full = [int((math.pi/2)**n) for n in range(1, num_terms_4 + 1)]
s4 = sorted(list(set(s4_full))) 
max_l_4 = find_max_ap_length(s4)
print(f"Generated {len(s4)} unique values up to {s4[-1]}.")
print(f"Found maximum AP length: {max_l_4}")
print("CONCLUSION: Theory states the max AP length is finite. The property HOLDS.")
print("-" * 50)
