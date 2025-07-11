import sys

# If you are using python version < 3.8, the modular inverse function is not available for negative numbers.
# The modular inverse function is implemented as pow(x, -1, m) in Python 3.8+.
# For older versions, pow(x, m-2, m) for a prime modulus m would be used.
if sys.version_info < (3, 8):
    print("This script requires Python 3.8+ for the pow(x, -1, m) function.")
    sys.exit(1)

# Define the modulus
MOD = 65539

# Step 1: Calculate S mod MOD
# S = 2^32. MOD = 2^16 + 3.
# S = (2^16)^2. Since 2^16 = MOD - 3, we have 2^16 === -3 (mod MOD).
# So, S === (-3)^2 = 9 (mod MOD).
S_mod = pow(2, 32, MOD)

# Step 2: Calculate m mod MOD
# m = 4 * S * (S - 4)
m_mod = (4 * S_mod * (S_mod - 4)) % MOD

# Step 3: Calculate M mod MOD
# M = S * (S - 4)^2 * (S + 8) / 27
# This is equivalent to S * (S-4)^2 * (S+8) * (27^-1) mod MOD

# Calculate modular inverse of 27 mod MOD
# Since MOD is a prime number, we can use Fermat's Little Theorem: a^(p-2) === a^-1 (mod p)
# Or in Python 3.8+, pow(27, -1, MOD)
inv_27 = pow(27, -1, MOD)

# Calculate numerator of M mod MOD
M_numerator_mod = (S_mod * pow(S_mod - 4, 2, MOD) * (S_mod + 8)) % MOD

# Calculate M mod MOD
M_mod = (M_numerator_mod * inv_27) % MOD

# Step 4: Calculate (m + M) mod MOD
result = (m_mod + M_mod) % MOD

# Print the numbers in the final equation
print(f"S mod {MOD} = {S_mod}")
print(f"m mod {MOD} = {m_mod}")
print(f"M mod {MOD} = {M_mod}")
print(f"(m + M) mod {MOD} = {result}")
