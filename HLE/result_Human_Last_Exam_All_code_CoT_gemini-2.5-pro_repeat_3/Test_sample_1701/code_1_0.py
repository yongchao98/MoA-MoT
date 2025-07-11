from fractions import Fraction

# Bloom filter parameters
m = 4  # number of bits in B
k = 2  # number of hash functions
n = 2  # number of inputs

# 1. Calculate p0: the probability a bit in B is 0.
# A bit is 0 if no hash function for any input maps to it.
# The formula is (1 - 1/m)^(k*n)
p0 = Fraction(1 - Fraction(1, m)) ** (k * n)

# 2. Calculate p1: the probability a bit in B is 1.
p1 = 1 - p0

# 3. Calculate the final probability for a bit in B' to be 1.
# B'[i] = B[j] XOR B[k]. This is 1 if B[j] != B[k].
# P(B'[i]=1) = P(B[j]=0, B[k]=1) + P(B[j]=1, B[k]=0)
# This simplifies to 2 * p0 * p1
final_prob = 2 * p0 * p1

print("This script calculates the probability that a bit in B' is 1.")
print(f"Given parameters: m={m}, k={k}, n={n}")
print("-" * 30)

# Print the equation for p0
print(f"Step 1: Calculate p0 (probability a bit in B is 0)")
print(f"p0 = (1 - 1/{m}) ** ({k} * {n})")
print(f"p0 = ({1-Fraction(1,m)}) ** {k*n}")
print(f"p0 = {p0.numerator}/{p0.denominator}")
print("")

# Print the equation for p1
print(f"Step 2: Calculate p1 (probability a bit in B is 1)")
print(f"p1 = 1 - p0 = 1 - {p0.numerator}/{p0.denominator}")
print(f"p1 = {p1.numerator}/{p1.denominator}")
print("")

# Print the final equation
print(f"Step 3: Calculate P(B'[i] = 1)")
print(f"P(B'[i] = 1) = 2 * p0 * p1")
print(f"P(B'[i] = 1) = 2 * ({p0.numerator}/{p0.denominator}) * ({p1.numerator}/{p1.denominator})")
print(f"P(B'[i] = 1) = {final_prob.numerator}/{final_prob.denominator}")
print("-" * 30)

# Print the final answer
print(f"The final probability is {float(final_prob)}")