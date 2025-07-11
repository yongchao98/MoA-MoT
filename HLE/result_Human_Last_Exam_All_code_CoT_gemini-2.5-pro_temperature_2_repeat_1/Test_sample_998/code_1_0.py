import math
from collections import Counter

def get_prime_factorization(n):
    """Computes the prime factorization of a number n."""
    factors = Counter()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    return factors

# Calculate the components of the sum S, as derived in the plan.
# S = 9 * S_b * S_eg

# S_b is the sum over the exponents of the prime 29.
# The summation becomes sum_{j=1 to 60} j.
S_b = sum(range(1, 60 + 1))

# S_eg is the sum over the exponents of primes 59 and 79.
# It splits based on the parity of the exponents.
sum_e_even = sum(80 - e for e in range(0, 80, 2))
sum_e_odd = sum(80 - e for e in range(1, 80, 2))
sum_g_even = sum(30 - g for g in range(0, 30, 2))
sum_g_odd = sum(30 - g for g in range(1, 30, 2))

S_eg = sum_e_even * sum_g_even + sum_e_odd * sum_g_odd

# The total sum S includes a factor of 9 from the tau(2^8) term.
S = 9 * S_b * S_eg

print(f"Let S be the sum of f(k) across all divisors of the given number.")
print("The value of S can be calculated as S = 9 * S_b * S_eg, where:")
print(f"S_b = {S_b}")
print(f"S_eg = ({sum_e_even} * {sum_g_even}) + ({sum_e_odd} * {sum_g_odd}) = {S_eg}")
print(f"S = 9 * {S_b} * {S_eg} = {S}")
print("-" * 30)

# To find the number of divisors of S, we find its prime factorization.
print("Prime factorization of S components:")
factors_9 = get_prime_factorization(9)
factors_sb = get_prime_factorization(S_b)
factors_seg = get_prime_factorization(S_eg)
print(f"  9 = {' * '.join([f'{p}^{e}' for p, e in sorted(factors_9.items())])}")
print(f"  {S_b} = {' * '.join([f'{p}^{e}' for p, e in sorted(factors_sb.items())])}")
print(f"  {S_eg} = {' * '.join([f'{p}^{e}' for p, e in sorted(factors_seg.items())])}")

# Combine the prime factors to get the factorization of S.
total_factors = factors_9 + factors_sb + factors_seg

S_factorization_str = " * ".join([f"{p}^{e}" for p, e in sorted(total_factors.items())])
print(f"\nThe prime factorization of S is: {S_factorization_str}")
print("-" * 30)

# Calculate the number of divisors using the exponents.
print("The number of divisors of S is calculated from the exponents of its prime factors:")
tau_eq_parts = []
num_divisors = 1
for p, e in sorted(total_factors.items()):
    tau_eq_parts.append(f"({e}+1)")
    num_divisors *= (e + 1)

final_equation = " * ".join(tau_eq_parts)
print(f"Number of divisors = {final_equation}")
print(f"Result = {num_divisors}")
