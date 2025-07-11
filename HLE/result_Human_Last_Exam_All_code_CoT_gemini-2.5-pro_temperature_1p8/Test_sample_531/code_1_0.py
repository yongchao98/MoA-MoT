# The modulus p
p = 65539

# S is 2^32. We compute S mod p.
# S = 2^32 = (2^16)^2 = 65536^2
# 65536 = 65539 - 3, so 65536 is congruent to -3 mod 65539.
# S is congruent to (-3)^2 = 9 mod 65539.
S_mod_p = pow(2, 32, p)

# Calculate m mod p
# m = 4*S*(S-4)
m_mod_p = (4 * S_mod_p * (S_mod_p - 4)) % p

# Calculate M mod p
# M = S * (S+8)/3 * ((S-4)/3)^2 = S * (S+8) * (S-4)^2 / 9
# To compute this modulo p, we can write M*9 = S*(S+8)*(S-4)^2
# (M mod p) * 9 = (S_mod_p * (S_mod_p+8) * (S_mod_p-4)^2) mod p
# We solve for M mod p by multiplying by the modular inverse of 9.
numerator_M_mod_p = (S_mod_p * (S_mod_p + 8) * pow(S_mod_p - 4, 2, p)) % p
inv_9_mod_p = pow(9, -1, p)
M_mod_p = (numerator_M_mod_p * inv_9_mod_p) % p

# The final result is (m + M) mod p.
result = (m_mod_p + M_mod_p) % p

# Output the numbers in the final equation as requested
# It will print "(180 + 425) % 65539 = 605"
print(f"({m_mod_p} + {M_mod_p}) % {p} = {result}")
