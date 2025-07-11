import sys

# It's recommended to run this code with python3.
# Let's define the modulus P
P = 65539

# Let's define S = 2^32. We will compute its value modulo P.
# S = 2^32 = (2^16)^2 = (65536)^2.
# Since 65539 is the modulus, 65536 = 65539 - 3.
# So, S_mod_P = (P-3)^2 mod P = (-3)^2 mod P = 9.
S_mod_P = pow(2, 32, P)
print(f"Let S = 2^32 and P = 65539.")
print(f"The value of S mod P is: {S_mod_P}")
print("-" * 30)

# Calculate m mod P
# m = S * (8*S - 48)
# We use modular arithmetic for the calculation.
m_mod_P = (S_mod_P * (8 * S_mod_P - 48)) % P
print("The second smallest value, m, is given by the formula m = S * (8*S - 48).")
print(f"m mod P = (S_mod_P * (8*S_mod_P - 48)) mod P")
print(f"m mod P = ({S_mod_P} * (8*{S_mod_P} - 48)) mod P")
print(f"m mod P = ({S_mod_P} * ({8*S_mod_P} - 48)) mod P")
print(f"m mod P = ({S_mod_P} * ({8*S_mod_P - 48})) mod P")
print(f"m mod P = {m_mod_P}")
print("-" * 30)

# Calculate M mod P
# M = S * (S-4)^2 * (S+8) / 27
# We compute (numerator * modular_inverse_of_denominator) mod P.
numerator = (S_mod_P * pow(S_mod_P - 4, 2, P) * (S_mod_P + 8)) % P

# To compute M mod P, we need the modular inverse of 27.
# Python's pow(base, exp, mod) can handle negative exponents for modular inverse.
inv_27 = pow(27, -1, P)
M_mod_P = (numerator * inv_27) % P

print("The second largest value, M, is given by the formula M = S * (S-4)^2 * (S+8) / 27.")
print("M mod P = (S_mod_P * (S_mod_P - 4)^2 * (S_mod_P + 8) * (27^-1)) mod P")
print(f"S_mod_P - 4 = {S_mod_P-4}")
print(f"S_mod_P + 8 = {S_mod_P+8}")
print(f"27^-1 mod P = {inv_27}")
print(f"M mod P = ({S_mod_P} * {(S_mod_P - 4)**2} * {S_mod_P + 8} * {inv_27}) mod P")
print(f"M mod P = {M_mod_P}")
print("-" * 30)


# Final result
result = (m_mod_P + M_mod_P) % P
print("The final result is (m + M) mod P.")
print(f"(m + M) mod P = ({m_mod_P} + {M_mod_P}) mod P = {result}")
