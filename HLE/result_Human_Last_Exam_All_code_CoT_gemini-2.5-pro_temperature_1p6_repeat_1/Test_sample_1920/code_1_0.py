import math

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
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

def phi(n):
    """Calculates Euler's totient function phi(n)."""
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

# Main part of the script
d = 53599
k = 6

print(f"Finding the number of primitive Dirichlet characters of conductor d={d} and order k={k}.")

# Step 1: Factorize d and verify conditions
factors = get_prime_factorization(d)
prime_factors = list(factors.keys())
print(f"The conductor d is {d}, with prime factorization: {' * '.join(map(str, prime_factors))}.")
num_factors = len(prime_factors)
print(f"Since d is square-free, a character is primitive if and only if all its {num_factors} component characters are non-principal (order > 1).")

# Step 2: Count available characters for each component
print(f"The order of the character is the lcm of component orders. For the lcm to be {k}, component orders must divide {k}.")
orders_to_consider = [o for o in [2, 3, 6] if k % o == 0] # Orders {2,3,6} are divisors of k=6
num_chars_by_order = {o: phi(o) for o in orders_to_consider}

print("\nFor each prime factor p, the number of primitive characters of order m > 1 is phi(m), provided m divides p-1.")
n2 = num_chars_by_order.get(2, 0)
n3 = num_chars_by_order.get(3, 0)
n6 = num_chars_by_order.get(6, 0)
print(f"Number of choices for a component character of order 2: phi(2) = {n2}")
print(f"Number of choices for a component character of order 3: phi(3) = {n3}")
print(f"Number of choices for a component character of order 6: phi(6) = {n6}")

# Step 3: Calculate the final number
total_choices_per_factor = n2 + n3 + n6
print(f"\nTotal choices for each component (orders dividing {k} but > 1): {total_choices_per_factor}")

total_combinations = total_choices_per_factor**num_factors
sub_lcm_2 = n2**num_factors
sub_lcm_3 = n3**num_factors

result = total_combinations - sub_lcm_2 - sub_lcm_3

print("\nUsing inclusion-exclusion, the total number is calculated as:")
print(f"({total_choices_per_factor})^{num_factors} - {n2}^{num_factors} - {n3}^{num_factors} = {result}")
