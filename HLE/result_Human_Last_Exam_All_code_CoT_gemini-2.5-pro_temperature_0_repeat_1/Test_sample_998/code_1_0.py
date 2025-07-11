import math

def get_prime_factorization_exponents(n):
    """
    Calculates the prime factorization of n and returns a list of exponents
    for prime factors sorted in ascending order.
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    
    # Return exponents of prime factors sorted by the prime base
    return [factors[p] for p in sorted(factors.keys())]

# Calculate the sum S based on the derived formula.
# Let N = 2^8 * 29^59 * 59^79 * 79^29.
# S = sum_{d|N} f(d) = sum_{k|N, k=1(mod 4)} tau(N/k).
# k = 29^b * 59^c * 79^e, where c+e is even.
# tau(N/k) = 9 * (60-b) * (80-c) * (30-e).

# Sum over b from 0 to 59
sum_b = sum(60 - b for b in range(60))

# Sums over c from 0 to 79
sum_c_even = sum(80 - c for c in range(0, 80, 2))
sum_c_odd = sum(80 - c for c in range(1, 80, 2))

# Sums over e from 0 to 29
sum_e_even = sum(30 - e for e in range(0, 30, 2))
sum_e_odd = sum(30 - e for e in range(1, 30, 2))

# The sum over c and e
sum_ce = sum_c_even * sum_e_even + sum_c_odd * sum_e_odd

# The total sum S
S = 9 * sum_b * sum_ce

# Get the exponents of the prime factorization of S
exponents = get_prime_factorization_exponents(S)

# Calculate the number of divisors of S
num_divisors = 1
equation_parts = []
for exp in exponents:
    term = exp + 1
    num_divisors *= term
    equation_parts.append(str(term))

# Format the output string as requested
equation_str = " * ".join(equation_parts)
print(f"The number of divisors is the result of the following product:")
print(f"{equation_str} = {num_divisors}")