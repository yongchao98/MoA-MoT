import math

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
    factors = {}
    d = 2
    temp_n = n
    # Handle factor of 2
    while temp_n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        temp_n //= 2
    # Handle odd factors
    d = 3
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 2
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

# Announce the goal of the script
print("This script calculates the highest possible order for the inertial quotient E.")
print("=" * 75)

# Step 1: Theoretical foundation
print("Step 1: Set up the problem based on block theory.")
p = 2
order_of_defect_group = 16
n = int(math.log(order_of_defect_group, p))
q = p

print(f"The defect group D is elementary abelian of order {order_of_defect_group} = {p}^{n}.")
print(f"This means D is a {n}-dimensional vector space over the field F_{p}.")
print(f"The automorphism group, Aut(D), is the general linear group GL({n}, {p}).")
print(f"The inertial quotient E is a p'-group (order is odd) and a subgroup of Aut(D).")
print(f"The highest possible order for E is the largest odd divisor of |GL({n}, {p})|.")
print("-" * 75)

# Step 2: Calculate the order of GL(n, q)
print(f"Step 2: Calculate the order of GL({n}, {p}).")
print("The formula for |GL(n, q)| is: (q^n - 1) * (q^n - q) * ... * (q^n - q^(n-1)).")

order_gl = 1
term_values = []
math_terms = []
for i in range(n):
    term = q**n - q**i
    order_gl *= term
    term_values.append(str(term))
    math_terms.append(f"({q}^{n} - {q}^{i})")

print(f"|GL({n}, {q})| = {' * '.join(math_terms)}")
print(f"           = ({16} - {1}) * ({16} - {2}) * ({16} - {4}) * ({16} - {8})")
print(f"           = {' * '.join(term_values)}")
print(f"           = {order_gl}")
print("-" * 75)

# Step 3: Find the largest odd divisor
print(f"Step 3: Find the largest odd divisor of {order_gl}.")
factors = get_prime_factorization(order_gl)

# Format the prime factorization string
factor_str_parts = []
for base, exp in sorted(factors.items()):
    factor_str_parts.append(f"{base}^{exp}" if exp > 1 else str(base))
print(f"The prime factorization of {order_gl} is: {' * '.join(factor_str_parts)}")

# Calculate the odd part and format the explanation string
odd_part = 1
odd_factors_symbolic = []
odd_factors_numeric = []
for base, exp in sorted(factors.items()):
    if base != p:
        odd_part *= (base ** exp)
        odd_factors_symbolic.append(f"{base}^{exp}" if exp > 1 else str(base))
        odd_factors_numeric.append(str(base**exp))

print("The largest odd divisor is the product of all prime factors other than 2.")
print(f"Highest order for E = {' * '.join(odd_factors_symbolic)}")
if len(odd_factors_numeric) > 1:
    print(f"                     = {' * '.join(odd_factors_numeric)}")
print(f"                     = {odd_part}")
print("=" * 75)
print(f"The highest order that E can have is {odd_part}.")
