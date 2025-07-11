import collections

def get_prime_factorization(n):
    """
    Returns a dictionary of the prime factorization of n.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

# --- Step 1: Define exponents from the problem statement ---
a2 = 8
a29 = 59
a59 = 79
a79 = 29

# --- Step 2: Calculate the components of the sum S ---

# Component for prime 29 (exponent b)
# This is the sum of integers from 1 to (a29 + 1)
sum_b = (a29 + 1) * (a29 + 2) // 2

# Components for prime 59 (exponent c)
sum_c_even = sum(a59 + 1 - c for c in range(0, a59 + 1, 2))
sum_c_odd = sum(a59 + 1 - c for c in range(1, a59 + 1, 2))

# Components for prime 79 (exponent e)
sum_e_even = sum(a79 + 1 - e for e in range(0, a79 + 1, 2))
sum_e_odd = sum(a79 + 1 - e for e in range(1, a79 + 1, 2))

# Combined component for c and e
# This sum is over c,e with the same parity (c+e is even)
sum_ce = sum_c_even * sum_e_even + sum_c_odd * sum_e_odd

# The leading factor from tau(N/k) related to the prime 2
factor_9 = a2 + 1

# Calculate the final value of S
S_val = factor_9 * sum_b * sum_ce

# --- Step 3: Find the prime factorization of S ---
total_factors = collections.defaultdict(int)

# Combine prime factors from all calculated components
for p, exp in get_prime_factorization(factor_9).items():
    total_factors[p] += exp
for p, exp in get_prime_factorization(sum_b).items():
    total_factors[p] += exp
for p, exp in get_prime_factorization(sum_ce).items():
    total_factors[p] += exp

# --- Step 4: Calculate the number of divisors of S ---
num_divisors = 1
explanation_parts = []
# Sort by prime for a consistent output format
sorted_primes = sorted(total_factors.keys())
for p in sorted_primes:
    exp = total_factors[p]
    num_divisors *= (exp + 1)
    explanation_parts.append(f"({exp} + 1)")

# --- Step 5: Print the results ---
print(f"The calculated value of S is the product of three main components:")
print(f"S = {factor_9} * {sum_b} * {sum_ce} = {S_val}")
print("-" * 20)
factor_string = " * ".join([f"{p}^{exp}" for p, exp in sorted(total_factors.items())])
print(f"The prime factorization of S is: {factor_string}")
print("-" * 20)
print(f"The number of divisors of S is the product of (exponent + 1) for each prime factor:")
equation = " * ".join(explanation_parts)
print(f"Number of divisors = {equation}")
print(f"Total number of divisors = {num_divisors}")
