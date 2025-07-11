import math
from functools import reduce

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors of n.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def get_all_factors(prime_factors):
    """
    Returns all factors from a dictionary of prime factors.
    """
    primes = list(prime_factors.keys())
    exponents = list(prime_factors.values())
    
    factors = [1]
    for i in range(len(primes)):
        new_factors = []
        for j in range(1, exponents[i] + 1):
            for factor in factors:
                new_factors.append(factor * (primes[i] ** j))
        factors.extend(new_factors)
    return sorted(factors)

# Step 1: Define the initial set and find its product's prime factorization
initial_set = list(range(1, 17))
primes_in_16_factorial = {}
for i in initial_set:
    factors = get_prime_factorization(i)
    for p, e in factors.items():
        primes_in_16_factorial[p] = primes_in_16_factorial.get(p, 0) + e

# Step 2: Determine the required exponents for the replacement product
replacement_product_exponents = {}
print("Prime factorization of 16! (exponents):")
for p, e in sorted(primes_in_16_factorial.items()):
    print(f"  Prime {p}: exponent {e}")
    if e % 4 != 0:
        replacement_product_exponents[p] = (4 - (e % 4)) % 4

print("\nRequired exponents (mod 4) for the product of replacement numbers k'1*k'2:")
# We are replacing 11 and 13.
# The exponents of 11 and 13 in 16! are 1.
# So we need to adjust the required exponents for these primes.
replacement_product_exponents[11] = (replacement_product_exponents.get(11, 0) + 1) % 4
replacement_product_exponents[13] = (replacement_product_exponents.get(13, 0) + 1) % 4
for p, e in sorted(replacement_product_exponents.items()):
    if e > 0:
        print(f"  For prime {p}, total exponent must be {e} (mod 4)")

# Step 3: Calculate the smallest integer product for the replacement pair
replacement_product = 1
for p, e in replacement_product_exponents.items():
    if e > 0:
        replacement_product *= (p ** e)

print(f"\nSmallest product of replacement numbers k'1 * k'2 is: {replacement_product}")

# Step 4: Find the factor pair of this product that minimizes the maximum value
factors_of_product = get_all_factors(get_prime_factorization(replacement_product))
best_pair = (1, replacement_product)
min_max_val = replacement_product
sqrt_prod = math.isqrt(replacement_product)

# Find factors around the square root
for factor in reversed(factors_of_product):
    if factor < sqrt_prod:
        break
    if replacement_product % factor == 0:
        other_factor = replacement_product // factor
        
        # Check if factors are in the remaining set
        remaining_set = set(range(1, 17)) - {11, 13}
        if factor in remaining_set or other_factor in remaining_set:
            continue
            
        current_max = max(factor, other_factor)
        if current_max < min_max_val:
            min_max_val = current_max
            best_pair = (other_factor, factor)

k1_prime, k2_prime = best_pair
final_N = max(k1_prime, k2_prime)

final_set = sorted(list((set(initial_set) - {11, 13}) | {k1_prime, k2_prime}))

print(f"\nOptimal replacement pair (k'1, k'2) is: {best_pair}")
print(f"The set of 16 distinct integers is: {final_set}")
print(f"The maximum value in this set is {max(final_set)}.")
print(f"The smallest N is {final_N}.")
