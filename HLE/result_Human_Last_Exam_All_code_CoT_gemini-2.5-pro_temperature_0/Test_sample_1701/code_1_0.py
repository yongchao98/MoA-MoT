from fractions import Fraction

# Parameters of the original Bloom filter B
m = 4  # number of bits
k = 2  # number of hash functions
n = 2  # number of inputs

# --- Step 1: Calculate the probability that a bit in B is 0 (p0) ---
# A bit in a Bloom filter remains 0 if no hash function maps any element to it.
# The probability that a single hash function for a single element does NOT map to a specific bit is (1 - 1/m).
# With k hash functions and n elements, there are k*n total hash computations.
# The probability that a bit is never set is (1 - 1/m)^(k*n).
p0 = Fraction(m - 1, m)**(k * n)

# --- Step 2: Calculate the probability that a bit in B is 1 (p1) ---
p1 = 1 - p0

# --- Step 3: Calculate the probability that a bit in the derived filter B' is 1 ---
# A bit B'[i] is 1 if B1[i] and B2[i] are different (B1[i] XOR B2[i] = 1).
# P(B'[i]=1) = P(B1[i]=0 and B2[i]=1) + P(B1[i]=1 and B2[i]=0)
# Assuming the states of bits in B are independent, this is:
# P(B'[i]=1) = p0 * p1 + p1 * p0 = 2 * p0 * p1
final_prob = 2 * p0 * p1

# --- Print the results with the full calculation ---
print("This script calculates the probability that a randomly chosen bit in B' is 1.\n")
print(f"Given parameters: m={m}, k={k}, n={n}\n")

print("--- Calculation Steps ---")

print("1. First, we find the probability (p0) that a bit in the original filter B is 0.")
print(f"   p0 = (1 - 1/m)^(k*n) = (1 - 1/{m})^({k}*{n})")
print(f"   p0 = ({p0.numerator}/{p0.denominator})\n")

print("2. Next, we find the probability (p1) that a bit in B is 1.")
print(f"   p1 = 1 - p0 = 1 - {p0.numerator}/{p0.denominator}")
print(f"   p1 = {p1.numerator}/{p1.denominator}\n")

print("3. Finally, we calculate the probability that a bit in B' is 1.")
print("   P(B'[i]=1) = P(B1[i] != B2[i]) = 2 * p0 * p1\n")

print("--- Final Equation ---")
print(f"P(B'[i]=1) = 2 * ({p0.numerator}/{p0.denominator}) * ({p1.numerator}/{p1.denominator})")
print(f"P(B'[i]=1) = {final_prob.numerator} / {final_prob.denominator}")
print(f"As a decimal, this is approximately: {float(final_prob)}")