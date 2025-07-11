import math

def prime_factorize(n):
    """
    Returns a list of unique prime factors of n.
    """
    factors = []
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.append(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.append(temp)
    return factors

def gcd(a, b):
    """
    Computes the greatest common divisor of a and b.
    """
    while b:
        a, b = b, a % b
    return a

def product(lst):
    """
    Computes the product of elements in a list.
    """
    res = 1
    for x in lst:
        res *= x
    return res

# Given values for the problem
d = 53599
g = 6

# Step 1: Find prime factors of the conductor d
prime_factors = prime_factorize(d)

# Step 2: Use inclusion-exclusion to find the number of characters.
# The lcm of orders must be exactly 6.
# This means the lcm must be divisible by 2 and 3.

# Let N(k) be the number of tuples of non-trivial characters (c_p1, c_p2, ...)
# such that lcm(ord(c_p1), ord(c_p2), ...) divides k.
# N(k) = product over p of (gcd(k, p-1) - 1).

# Calculate N(6): Total combinations where lcm(orders) divides 6
n_div_6_list = [gcd(6, p - 1) - 1 for p in prime_factors]
N6 = product(n_div_6_list)

# Calculate N(3): Combinations where lcm(orders) divides 3 (lacks factor of 2)
n_div_3_list = [gcd(3, p - 1) - 1 for p in prime_factors]
N3 = product(n_div_3_list)

# Calculate N(2): Combinations where lcm(orders) divides 2 (lacks factor of 3)
n_div_2_list = [gcd(2, p - 1) - 1 for p in prime_factors]
N2 = product(n_div_2_list)

# The number of characters with order lcm being exactly 6 is N6 - N3 - N2
# The intersection N(1) is 0 because characters are non-trivial (order > 1).
result = N6 - N3 - N2

# Step 3: Print the final equation with all numbers.
print(f"The conductor is d = {d}, which has prime factors: {prime_factors}.")
print(f"The required order is g = {g}.")
print("Let A be the number of primitive characters where the LCM of the orders of the component characters divides 6.")
print("Let B be the number where the LCM of orders divides 3.")
print("Let C be the number where the LCM of orders divides 2.")
print("By the principle of inclusion-exclusion, the desired number is N = A - B - C.")
print(f"A = {' * '.join(map(str, n_div_6_list))} = {N6}")
print(f"B = {' * '.join(map(str, n_div_3_list))} = {N3}")
print(f"C = {' * '.join(map(str, n_div_2_list))} = {N2}")
print(f"\nThe number of such characters is N = {N6} - {N3} - {N2} = {result}")
