import math

# Parameters from the IEEE 754 64-bit example
m = 52
e = 11

# --- Case 1: Bias as defined in the IEEE 754 example ---
# B = 2^(e-1) - 1
B1 = 2**(e - 1) - 1
# Total number of representable values in [0, 1]
N1 = B1 * (2**m) + 1
# Minimum bits required is the ceiling of log2(N)
bits1 = math.ceil(math.log2(N1))

print(f"Using the example bias B = 2^(e-1) - 1 = {B1}:")
print(f"Number of values N = B * 2^m + 1 = {B1} * 2^{m} + 1")
print(f"log2(N) = {math.log2(N1)}")
print(f"Fewest bits required = ceil(log2(N)) = {bits1}")
print(f"This equals m + e - 1 = {m} + {e} - 1 = {m + e - 1}\n")


# --- Case 2: Bias based on an alternative common standard ---
# B = 2^(e-1)
B2 = 2**(e - 1)
# Total number of representable values in [0, 1]
N2 = B2 * (2**m) + 1
# Minimum bits required
bits2 = math.ceil(math.log2(N2))

print(f"Using an alternative bias B = 2^(e-1) = {B2}:")
print(f"Number of values N = B * 2^m + 1 = {B2} * 2^{m} + 1")
print(f"log2(N) = {math.log2(N2)}")
print(f"Fewest bits required = ceil(log2(N)) = {bits2}")
print(f"This equals m + e = {m} + {e} = {m + e}")
