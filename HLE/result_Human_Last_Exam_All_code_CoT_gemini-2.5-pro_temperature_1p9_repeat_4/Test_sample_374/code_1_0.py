import math

def get_prime_factorization(num):
    """Returns the prime factorization of num as a dictionary."""
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

# 1. Define the parameters from the problem statement
n = 4 # Dimension of the vector space, from D being order 2^4
p = 2 # Characteristic of the field k

print("The inertial quotient E is a p'-subgroup of Out(D).")
print(f"Here, p=2, so E has odd order. D is elementary abelian of order 16=2^{n}, so Out(D) is GL({n},F_{p}).")
print("We need to find the odd part of the order of GL(4, 2).\n")

# 2. Calculate the order of GL(n, p)
order = 1
# Formula is product from i=0 to n-1 of (p^n - p^i)
print(f"The order of GL({n}, {p}) is calculated as:")
calc_str_list = []
val_str_list = []
for i in range(n):
    term = (p**n - p**i)
    order *= term
    calc_str_list.append(f"({p}^{n} - {p}^{i})")
    val_str_list.append(str(term))

print(f"|GL({n}, {p})| = " + " * ".join(calc_str_list))
print(f"             = " + " * ".join(val_str_list))
print(f"             = {order}\n")

# 3. Find the prime factorization of the order
prime_factors = get_prime_factorization(order)
print(f"The prime factorization of {order} is:")
factor_strings = []
for prime, exp in sorted(prime_factors.items()):
    if exp > 1:
        factor_strings.append(f"{prime}^{exp}")
    else:
        factor_strings.append(f"{prime}")
print(f"{order} = {' * '.join(factor_strings)}\n")


# 4. Calculate the largest odd divisor (the p'-part)
largest_odd_divisor = 1
odd_factor_strings_exp = []
odd_factor_strings_val = []

for prime, exp in prime_factors.items():
    if prime != p: # In this case, p=2
        term_val = prime**exp
        largest_odd_divisor *= term_val
        if exp > 1:
            odd_factor_strings_exp.append(f"{prime}^{exp}")
        else:
            odd_factor_strings_exp.append(f"{prime}")
        odd_factor_strings_val.append(str(term_val))

print("The highest order for E is the odd part of this number.")
print("Final Calculation:")
print(f"{' * '.join(odd_factor_strings_exp)} = {' * '.join(odd_factor_strings_val)} = {largest_odd_divisor}")