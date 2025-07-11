from fractions import Fraction

# Parameters of the Bloom filter
m = 4  # number of bits in B
k = 2  # number of hash functions
n = 2  # number of inputs

# Total number of hash computations
H = k * n

# The probability that a randomly chosen bit in B' is 1 is the same as
# P(b_i XOR b_j = 1) for any two distinct bits b_i, b_j from B.
# This can be calculated with the formula:
# P(B'[i]=1) = 2 * (P(b_i = 0) - P(b_i = 0, b_j = 0))

# --- Step 1: Calculate P(b_i = 0) ---
# This is the probability that a single bit in B is 0, which means all H
# hashes missed this bit. The probability for a single hash to miss is (m-1)/m.
p_bi_is_0 = Fraction(m - 1, m) ** H

# --- Step 2: Calculate P(b_i = 0, b_j = 0) ---
# This is the probability that two distinct bits in B are both 0. This means all H
# hashes missed both bits. The probability for a single hash to miss both is (m-2)/m.
p_bi_and_bj_are_0 = Fraction(m - 2, m) ** H

# --- Step 3: Calculate the final probability ---
prob = 2 * (p_bi_is_0 - p_bi_and_bj_are_0)

# --- Print the step-by-step calculation ---
print("This script calculates the probability that a randomly chosen bit in B' is 1.")
print("-" * 60)
print(f"Parameters: m={m} (bits), k={k} (hash functions), n={n} (inputs)")
print(f"Total hash computations H = k * n = {H}")
print("\nFinal Equation: P(B'[i]=1) = 2 * (P(b_i = 0) - P(b_i = 0, b_j = 0))")
print("")

print("Calculation Steps:")
print(f"1. P(b_i = 0) = (({m-1})/{m})^{H} = ({p_bi_is_0.numerator})/({p_bi_is_0.denominator})")
print(f"2. P(b_i = 0, b_j = 0) = (({m-2})/{m})^{H} = ({p_bi_and_bj_are_0.numerator})/({p_bi_and_bj_are_0.denominator})")
print("")

# To show the subtraction with a common denominator
common_denominator = p_bi_is_0.denominator
val1_num = p_bi_is_0.numerator
val2_num = p_bi_and_bj_are_0.numerator * (common_denominator // p_bi_and_bj_are_0.denominator)
sub_result_num = val1_num - val2_num

print("3. Substitute these values into the formula:")
print(f"   P(B'[i]=1) = 2 * ( {p_bi_is_0} - {p_bi_and_bj_are_0} )")
print(f"              = 2 * ( {val1_num}/{common_denominator} - {val2_num}/{common_denominator} )")
print(f"              = 2 * ( {sub_result_num}/{common_denominator} )")
print(f"              = {2 * sub_result_num}/{common_denominator}")
print(f"              = {prob.numerator}/{prob.denominator}")
print("")
print("Final Answer:")
print(f"The probability that a randomly chosen bit in B' is 1 is {prob.numerator}/{prob.denominator}, or approximately {float(prob):.7f}.")
