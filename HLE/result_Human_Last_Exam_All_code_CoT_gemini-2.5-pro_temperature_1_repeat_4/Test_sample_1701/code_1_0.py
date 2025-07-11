import math

# Parameters of the Bloom filter as defined in the problem
m = 4  # Total bits in the original filter B
k = 2  # Number of hash functions
n = 2  # Number of inputs

# Total number of hash function applications
N = k * n

# The problem asks for P(B'[i] = 1), where B'[i] = B_1[i] XOR B_2[i].
# This is equivalent to finding the probability that two corresponding bits
# from the two halves of B are different. By symmetry, this is the same
# as finding the probability that any two distinct bits in B, say b_j and b_l,
# have different values.

# The probability that b_j != b_l is P(b_j=1, b_l=0) + P(b_j=0, b_l=1).
# By symmetry, this is 2 * P(b_j=1, b_l=0).
# We can find P(b_j=1, b_l=0) from P(b_l=0) - P(b_j=0, b_l=0).

# Step 1: Calculate the probability that a single specific bit b_l is 0.
# This occurs if none of the N hash applications fall on position l.
# For one hash, P(not hitting l) = (m-1)/m. For N hashes, it's ((m-1)/m)^N.
prob_bl_is_0_num = (m - 1)**N
prob_bl_is_0_den = m**N

# Step 2: Calculate the probability that two specific bits b_j and b_l are both 0.
# This occurs if none of the N hash applications fall on position j or l.
# For one hash, P(not hitting j or l) = (m-2)/m. For N hashes, it's ((m-2)/m)^N.
prob_bj0_and_bl0_num = (m - 2)**N
prob_bj0_and_bl0_den = m**N

# Step 3: Calculate P(b_j=1, b_l=0).
# P(b_j=1, b_l=0) = P(b_l=0) - P(b_j=0, b_l=0)
# To subtract the fractions, we need a common denominator.
common_den = prob_bl_is_0_den
# Adjust the numerator of the second term
adjusted_bj0_bl0_num = prob_bj0_and_bl0_num * (common_den // prob_bj0_and_bl0_den)
prob_bj1_and_bl0_num = prob_bl_is_0_num - adjusted_bj0_bl0_num
prob_bj1_and_bl0_den = common_den

# Step 4: Calculate the final probability P(b_j XOR b_l = 1).
# P(XOR=1) = 2 * P(b_j=1, b_l=0)
prob_xor_1_num = 2 * prob_bj1_and_bl0_num
prob_xor_1_den = prob_bj1_and_bl0_den

# Step 5: Simplify the final fraction.
common_divisor = math.gcd(prob_xor_1_num, prob_xor_1_den)
final_num = prob_xor_1_num // common_divisor
final_den = prob_xor_1_den // common_divisor

# Print the step-by-step derivation with numbers
print("Derivation of the probability P(B'[i] = 1):")
print("-" * 40)
print(f"Let b_j and b_l be two distinct bits in the original filter B.")
print(f"P(B'[i] = 1) = P(b_j != b_l) = 2 * P(b_j=1, b_l=0)")
print(f"We use the formula: P(b_j=1, b_l=0) = P(b_l=0) - P(b_j=0, b_l=0)\n")

print(f"1. Probability that a bit b_l is 0:")
print(f"   P(b_l=0) = (({m}-1)/{m})^({k}*{n}) = ({m-1}/{m})^{N} = {prob_bl_is_0_num}/{prob_bl_is_0_den}")

print(f"\n2. Probability that two bits b_j and b_l are both 0:")
print(f"   P(b_j=0, b_l=0) = (({m}-2)/{m})^({k}*{n}) = ({m-2}/{m})^{N} = {prob_bj0_and_bl0_num}/{prob_bj0_and_bl0_den}")

print(f"\n3. Probability that b_j=1 and b_l=0:")
print(f"   P(b_j=1, b_l=0) = {prob_bl_is_0_num}/{prob_bl_is_0_den} - {adjusted_bj0_bl0_num}/{common_den} = {prob_bj1_and_bl0_num}/{prob_bj1_and_bl0_den}")

print(f"\n4. Final probability that the XOR is 1:")
print(f"   P(B'[i]=1) = 2 * ({prob_bj1_and_bl0_num}/{prob_bj1_and_bl0_den}) = {prob_xor_1_num}/{prob_xor_1_den}")

print(f"\nSimplified, the final probability is {final_num}/{final_den}")
print(f"As a decimal, the probability is {final_num/final_den}")
