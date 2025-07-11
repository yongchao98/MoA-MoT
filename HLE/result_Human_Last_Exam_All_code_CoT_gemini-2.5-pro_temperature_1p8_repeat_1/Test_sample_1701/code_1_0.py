import math
from fractions import Fraction

# --- Problem Parameters ---
m = 4  # Number of bits in the Bloom filter B
k = 2  # Number of hash functions
n = 2  # Number of inputs

# --- Step-by-step Calculation ---

print("Problem: Find the probability P(B'[i] = 1)")
print(f"where B' is derived from a Bloom filter B with m={m}, k={k}, n={n}.")
print("B'[i] = B_1[i] XOR B_2[i], where B_1 and B_2 are the two halves of B.")
print("We will calculate P(B'[1] = 1) = P(B[1] XOR B[3] = 1) = P(B[1] != B[3]).")
print("This probability is P(B[1]=0, B[3]=1) + P(B[1]=1, B[3]=0).\n")

# Total number of hashing operations
nk = n * k
print(f"1. Total number of hash operations = n * k = {n} * {k} = {nk}.\n")

# --- Probability of a single bit being 0 ---
# This occurs if none of the 'nk' hashes map to this bit.
# The probability for one hash to miss is (m-1)/m.
prob_one_hash_misses = Fraction(m - 1, m)
p_j_zero = prob_one_hash_misses ** nk

print("2. Calculate the probability that a single specific bit B[j] is 0.")
print(f"   P(B[j]=0) = ((m-1)/m)^(n*k)")
print(f"   P(B[j]=0) = (({m}-1)/{m})^({nk}) = ({p_j_zero.numerator}/{p_j_zero.denominator})\n")

# --- Probability of two specific bits being 0 ---
# This occurs if none of the 'nk' hashes map to either of these two bits.
# The probability for one hash to miss is (m-2)/m.
prob_one_hash_misses_two_bits = Fraction(m - 2, m)
p_jl_zero = prob_one_hash_misses_two_bits ** nk

print("3. Calculate the probability that two specific bits B[j] and B[l] are both 0.")
print(f"   P(B[j]=0, B[l]=0) = ((m-2)/m)^(n*k)")
print(f"   P(B[j]=0, B[l]=0) = (({m}-2)/{m})^({nk}) = ({p_jl_zero.numerator}/{p_jl_zero.denominator})\n")

# --- Probability of one bit being 1 and another being 0 ---
# P(B[1]=1, B[3]=0) = P(B[3]=0) - P(B[1]=0, B[3]=0)
# By symmetry, this is the same as P(B[1]=0, B[3]=1).
prob_one_one_three_zero = p_j_zero - p_jl_zero

print("4. Calculate the probability that B[1]=1 AND B[3]=0.")
print(f"   P(B[1]=1, B[3]=0) = P(B[3]=0) - P(B[1]=0, B[3]=0)")
print(f"   P(B[1]=1, B[3]=0) = {p_j_zero.numerator}/{p_j_zero.denominator} - {p_jl_zero.numerator}/{p_jl_zero.denominator} = {prob_one_one_three_zero.numerator}/{prob_one_one_three_zero.denominator}\n")

# --- Final Probability Calculation ---
# P(B'[1]=1) = P(B[1]=1, B[3]=0) + P(B[1]=0, B[3]=1)
# P(B'[1]=1) = 2 * P(B[1]=1, B[3]=0)
final_prob = 2 * prob_one_one_three_zero

print("5. Calculate the final probability P(B'[1]=1).")
print(f"   P(B'[1]=1) = P(B[1]=1, B[3]=0) + P(B[1]=0, B[3]=1)")
print(f"   The equation is: 2 * (P(B[j]=0) - P(B[j]=0, B[l]=0))")
print(f"   P(B'[1]=1) = 2 * ({p_j_zero.numerator}/{p_j_zero.denominator} - {p_jl_zero.numerator}/{p_jl_zero.denominator}) = 2 * ({prob_one_one_three_zero.numerator}/{prob_one_one_three_zero.denominator}) = {final_prob.numerator}/{final_prob.denominator}\n")

print("--- Final Answer ---")
print(f"The probability that a randomly chosen bit in B' is 1 is {final_prob.numerator}/{final_prob.denominator}")
print(f"As a decimal, this is: {float(final_prob)}")