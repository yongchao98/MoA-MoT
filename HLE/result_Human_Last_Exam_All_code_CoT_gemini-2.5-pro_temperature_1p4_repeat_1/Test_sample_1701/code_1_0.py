import fractions

# Step 1: Define problem parameters
n = 2  # number of inputs
k = 2  # number of hash functions
m = 4  # number of bits in the original Bloom filter B

# The problem asks for the probability P(B'[i] = 1) for a randomly chosen bit i.
# Due to symmetry, this is equal to P(B'[1] = 1).
# B'[1] = B[1] XOR B[3].
# P(B'[1] = 1) = P(B[1] != B[3]) = 2 * P(B[1]=0, B[3]=1)
# We will calculate P(B[1]=0, B[3]=1) = P(B[3]=1 | B[1]=0) * P(B[1]=0)

# Total number of hash operations
N = n * k

print(f"Solving for a Bloom filter with n={n}, k={k}, m={m}.")
print(f"The total number of hash operations is N = n * k = {N}.")
print("-" * 30)
print("The final probability is calculated with the formula:")
print(f"P(B'[i]=1) = 2 * (1 - (({m}-2)/({m}-1))^{N}) * (({m}-1)/{m})^{N}")
print("-" * 30)

# Step 2: Calculate P(B[1] = 0)
# This is the probability that none of the N hashes land on bit 1.
# The probability of one hash not landing on bit 1 is (m-1)/m.
p_b1_is_0_num = (m - 1)**N
p_b1_is_0_den = m**N
p_b1_is_0 = fractions.Fraction(p_b1_is_0_num, p_b1_is_0_den)

print("Calculation steps:")
print(f"1. Probability that bit B[1] is 0:")
print(f"   P(B[1]=0) = (({m}-1)/{m})^{N} = ({m-1}/{m})^{N} = {p_b1_is_0.numerator}/{p_b1_is_0.denominator}")

# Step 3: Calculate P(B[3]=1 | B[1]=0)
# Given B[1]=0, all N hashes landed in the other m-1 bits.
# This is 1 minus the probability that B[3] is also 0 in this new scenario.
# The probability for a hash to not hit bit 3 (out of m-1 options) is (m-2)/(m-1).
p_b3_is_0_cond_b1_is_0_num = (m - 2)**N
p_b3_is_0_cond_b1_is_0_den = (m - 1)**N
p_b3_is_1_cond_b1_is_0 = 1 - fractions.Fraction(p_b3_is_0_cond_b1_is_0_num, p_b3_is_0_cond_b1_is_0_den)

print(f"2. Conditional probability that B[3] is 1, given B[1] is 0:")
print(f"   P(B[3]=1 | B[1]=0) = 1 - (({m}-2)/({m}-1))^{N} = 1 - ({m-2}/{m-1})^{N} = {p_b3_is_1_cond_b1_is_0.numerator}/{p_b3_is_1_cond_b1_is_0.denominator}")

# Step 4: Calculate the joint probability P(B[1]=0, B[3]=1)
p_b1_0_and_b3_1 = p_b3_is_1_cond_b1_is_0 * p_b1_is_0
print(f"3. Joint probability of B[1]=0 and B[3]=1:")
print(f"   P(B[1]=0, B[3]=1) = {p_b3_is_1_cond_b1_is_0} * {p_b1_is_0} = {p_b1_0_and_b3_1.numerator}/{p_b1_0_and_b3_1.denominator}")

# Step 5: Calculate the final probability, P(B'[1] = 1)
final_prob = 2 * p_b1_0_and_b3_1
print(f"4. Final probability P(B'[i]=1):")
print(f"   P(B'[i]=1) = 2 * P(B[1]=0, B[3]=1) = 2 * {p_b1_0_and_b3_1} = {final_prob.numerator}/{final_prob.denominator}")

print("-" * 30)
print(f"The final probability that a randomly chosen bit in B' is 1 is {final_prob.numerator}/{final_prob.denominator}.")
<<<65/128>>>