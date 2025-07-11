import fractions

# Problem parameters
m = 4  # number of bits in the original Bloom filter B
k = 2  # number of hash functions
n = 2  # number of inputs

# The probability we want to find is P(B'[i] = 1).
# B'[i] is 1 if and only if B_1[i] and B_2[i] are different.
# Let b1 and b2 be the corresponding bits from the two halves of B.
# P(B'[i]=1) = P(b1=0, b2=1) + P(b1=1, b2=0)
# By symmetry, this is 2 * P(b1=0, b2=1).
# We can find P(b1=0, b2=1) using: P(b1=0, b2=1) = P(b1=0) - P(b1=0, b2=0).

print("This script calculates the probability that a randomly chosen bit in B' is 1.")
print("The final probability is calculated as 2 * (p0 - p00), where:")
print("p0: The probability that a single bit in the original Bloom filter is 0.")
print("p00: The probability that two specific bits in the original Bloom filter are both 0.")
print("-" * 50)

# Calculate k*n, the total number of hashings
kn = k * n

# Step 1: Calculate p0 = P(a specific bit is 0)
# This occurs if all kn hashings miss the bit.
# The probability a single hashing misses is (1 - 1/m).
# The probability all kn hashings miss is (1 - 1/m)^(kn).
p0_num = (m - 1)**kn
p0_den = m**kn
p0 = fractions.Fraction(p0_num, p0_den)

print("Step 1: Calculating p0")
print(f"p0 = (1 - 1/m)^(k*n) = (1 - 1/{m})^({k}*{n})")
print(f"p0 = (({m-1}/{m}))^{kn} = {(m-1)**kn}/{m**kn} = {p0}")
print("-" * 50)

# Step 2: Calculate p00 = P(two specific bits are both 0)
# This occurs if all kn hashings miss both bits.
# The probability a single hashing misses both is (1 - 2/m).
# The probability all kn hashings miss both is (1 - 2/m)^(kn).
p00_num = (m - 2)**kn
p00_den = m**kn
p00 = fractions.Fraction(p00_num, p00_den)

print("Step 2: Calculating p00")
print(f"p00 = (1 - 2/m)^(k*n) = (1 - 2/{m})^({k}*{n})")
print(f"p00 = (({m-2}/{m}))^{kn} = {(m-2)**kn}/{m**kn} = {p00}")
print("-" * 50)

# Step 3: Calculate the probability that one bit is 0 and the other is 1
p_01 = p0 - p00

print("Step 3: Calculating P(bit_1=0, bit_2=1)")
print(f"P(bit_1=0, bit_2=1) = p0 - p00 = {p0} - {p00} = {p_01}")
print("-" * 50)

# Step 4: Calculate the final probability P(XOR=1)
# P(XOR=1) = 2 * P(bit_1=0, bit_2=1)
final_prob = 2 * p_01

print("Step 4: Calculating the final probability P(B'[i]=1)")
print("P(B'[i]=1) = P(bit_1 != bit_2) = 2 * P(bit_1=0, bit_2=1)")
print(f"P(B'[i]=1) = 2 * {p_01} = {final_prob}")
print("-" * 50)

print(f"\nThe final probability that a randomly chosen bit in B' is 1 is {final_prob}.")
print(f"As a decimal, this is {float(final_prob)}.")
