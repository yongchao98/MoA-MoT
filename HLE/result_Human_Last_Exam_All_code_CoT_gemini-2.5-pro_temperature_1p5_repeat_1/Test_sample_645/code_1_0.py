import numpy as np
import math

def generate_primes(n):
    """Generate the first n primes."""
    primes = []
    num = 2
    if n > 0:
        primes.append(num)
    num = 3
    while len(primes) < n:
        is_prime = True
        # Check divisibility only up to the square root of num
        for p in primes:
            if p * p > num:
                break
            if num % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(num)
        num += 2 # Check only odd numbers
    return primes

def get_set_1(num_terms):
    """
    S = {sum_{k<=n} N_k} where N_k ~ Poi(1).
    Returns a realization of the set of unique values.
    """
    # Fix a seed for reproducibility
    np.random.seed(42)
    N_k = np.random.poisson(1, num_terms)
    s_n = np.cumsum(N_k)
    return sorted(list(set(s_n)))

def get_set_2(num_terms, k=4):
    """S = {n^k for n in N}"""
    return [n**k for n in range(1, num_terms + 1)]

def get_set_3(num_terms):
    """S = the set of primes"""
    return generate_primes(num_terms)

def get_set_4(num_terms):
    """S = {floor((pi/2)^n) for n in N}"""
    return [math.floor((math.pi/2)**n) for n in range(1, num_terms + 1)]

# --- Main execution ---
num_terms_display = 20

print("This script generates the first few elements of each set to illustrate their structure.")
print("The solution to the problem relies on the mathematical theory of Sidon sets.\n")

# Set 1
# Generate more raw terms to get a decent number of unique values
s1_raw_len = 30
s1 = get_set_1(s1_raw_len)
print(f"Set 1 (a realization of sums of Poisson variables, first {len(s1)} unique terms from {s1_raw_len} steps):")
print(s1)
if len(s1) > 1 and s1[-1] > 0:
    density = len(s1) / s1[-1]
    print(f"Structure: Dense. Approximate density: {len(s1)}/{s1[-1]} = {density:.2f}")
    print("Conclusion: Not a Sidon set. Has the property.\n")

# Set 2
s2 = get_set_2(num_terms_display)
print(f"Set 2 (k=4, first {num_terms_display} terms):")
print(s2)
if len(s2) > 1:
  ratios = [s2[i+1]/s2[i] for i in range(len(s2)-1)]
  print(f"Structure: Polynomial growth. Ratios of consecutive terms -> 1: {[f'{r:.2f}' for r in ratios[:5]]}...")
  print("Conclusion: Not a Sidon set. Has the property.\n")

# Set 3
s3 = get_set_3(num_terms_display)
print(f"Set 3 (primes, first {num_terms_display} terms):")
print(s3)
if len(s3) > 1:
  ratios = [s3[i+1]/s3[i] for i in range(len(s3)-1)]
  print(f"Structure: Sparse (density ~1/log N). Ratios -> 1: {[f'{r:.2f}' for r in ratios[:5]]}...")
  print("Conclusion: Is a Sidon set. Does not have the property.\n")

# Set 4
s4 = get_set_4(num_terms_display)
print(f"Set 4 (floor((pi/2)^n), first {num_terms_display} terms):")
print(s4)
if len(s4) > 1:
  # Filter out s_i=0 for ratio calculation, happens if (pi/2)^n < 1
  s4_positive = [s for s in s4 if s > 0]
  if len(s4_positive) > 1:
    ratios = [s4_positive[i+1]/s4_positive[i] for i in range(len(s4_positive)-1)]
    print(f"Structure: Lacunary (exponential growth). Ratios -> pi/2 ~ 1.57: {[f'{r:.2f}' for r in ratios[:5]]}...")
  print("Conclusion: Is a Sidon set. Does not have the property.\n")