from fractions import Fraction

# --- Problem Parameters ---
# m: number of bits in the Bloom filter B
# n: number of inputs inserted
# k: number of hash functions
m = 4
n = 2
k = 2

# Total number of hash function applications
nk = n * k

print("This script calculates the probability that a randomly chosen bit in B' is 1.")
print(f"Parameters: m={m} bits, n={n} inputs, k={k} hash functions.\n")

# --- Step 1: Calculate the probability that a single bit in B is 0 ---
# A specific bit remains 0 if none of the nk hash functions map to its position.
# The probability of a single hash not mapping to a bit is (m-1)/m.
# For nk independent hashes, this is ((m-1)/m)^nk.
p0_num = (m - 1)**nk
p0_den = m**nk
p0 = Fraction(p0_num, p0_den)

print("Step 1: Calculate P(bit=0), the probability a single bit in B is 0.")
print(f"P(bit=0) = (({m-1})/{m})^({n*k}) = {p0_num}/{p0_den}")
print(f"P(bit=0) = {p0}\n")


# --- Step 2: Calculate the probability that two specific bits in B are both 0 ---
# Two specific bits remain 0 if none of the nk hash functions map to either position.
# The probability of a single hash not mapping to either bit is (m-2)/m.
# For nk independent hashes, this is ((m-2)/m)^nk.
p00_num = (m - 2)**nk
p00_den = m**nk
p00 = Fraction(p00_num, p00_den)

print("Step 2: Calculate P(bit_j=0, bit_k=0), the probability two specific bits in B are 0.")
print(f"P(bit_j=0, bit_k=0) = (({m-2})/{m})^({n*k}) = {p00_num}/{p00_den}")
print(f"P(bit_j=0, bit_k=0) = {p00}\n")


# --- Step 3: Calculate the final probability for B'[i] = 1 ---
# P(B'[i]=1) = P(bit_j=1, bit_k=0) + P(bit_j=0, bit_k=1)
# This simplifies to 2 * (P(bit=0) - P(bit_j=0, bit_k=0))
prob_xor_one = 2 * (p0 - p00)

print("Step 3: Calculate the final probability P(B'[i]=1).")
print("The formula is: 2 * (P(bit=0) - P(bit_j=0, bit_k=0))")
print(f"P(B'[i]=1) = 2 * ({p0} - {p00})")
print(f"P(B'[i]=1) = 2 * ({p0 - p00})")
print(f"P(B'[i]=1) = {prob_xor_one}\n")

final_float = float(prob_xor_one)
print(f"The final probability is {prob_xor_one}, or approximately {final_float:.5f}.")