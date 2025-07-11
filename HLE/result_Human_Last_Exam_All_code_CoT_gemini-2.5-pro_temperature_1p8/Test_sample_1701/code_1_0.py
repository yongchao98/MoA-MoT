import fractions

# Step 1: Define the parameters of the Bloom filter problem.
m = 4  # Total number of bits in the filter B
k = 2  # Number of hash functions
n = 2  # Number of inputs inserted into the filter

# The derived array B' is formed by XORing pairs of bits from B.
# We want to find the probability that a bit in B' is 1.
# P(B'[i]=1) = P(B_x XOR B_y = 1), where B_x and B_y are distinct bits from B.
# This occurs when (B_x=1 AND B_y=0) or (B_x=0 AND B_y=1).
# So, P(B'[i]=1) = P(B_x=1, B_y=0) + P(B_x=0, B_y=1).

# Step 2: Calculate the total number of hash operations.
nk = n * k

print(f"The probability of a randomly chosen bit in B' being 1 is calculated as follows:\n")
print(f"P(B'[i]=1) = P(B_x=1, B_y=0) + P(B_x=0, B_y=1)\n")

# Step 3: Calculate the probability that two distinct bits, B_x and B_y, are both 0.
# This occurs if none of the nk hashes point to either position x or y.
# For a single hash, the probability of avoiding both is (m-2)/m.
prob_both_0_num = (m - 2)**nk
prob_both_0_den = m**nk
f_prob_both_0 = fractions.Fraction(prob_both_0_num, prob_both_0_den)

print(f"1. Probability that two specific bits (B_x, B_y) are both 0:")
print(f"   P(B_x=0, B_y=0) = (({m}-2)/{m})^({n}*{k}) = ({m-2}/{m})^{nk} = {f_prob_both_0.numerator}/{f_prob_both_0.denominator}\n")

# Step 4: Calculate the probability that a single bit, B_y, is 0.
# This occurs if none of the nk hashes point to position y.
# For a single hash, the probability of avoiding y is (m-1)/m.
prob_one_0_num = (m - 1)**nk
prob_one_0_den = m**nk
f_prob_one_0 = fractions.Fraction(prob_one_0_num, prob_one_0_den)

print(f"2. Probability that one specific bit (B_y) is 0:")
print(f"   P(B_y=0) = (({m}-1)/{m})^({n}*{k}) = ({m-1}/{m})^{nk} = {f_prob_one_0.numerator}/{f_prob_one_0.denominator}\n")

# Step 5: Calculate P(B_x=1, B_y=0) using the formula P(A and not B) = P(not B) - P(not A and not B).
# Let A be the event B_x=1 and B be the event B_y=1.
# P(B_x=1, B_y=0) = P(B_y=0) - P(B_x=0, B_y=0)
f_prob_one_1_one_0 = f_prob_one_0 - f_prob_both_0

print(f"3. Probability that B_x=1 and B_y=0:")
print(f"   P(B_x=1, B_y=0) = P(B_y=0) - P(B_x=0, B_y=0) = {f_prob_one_0.numerator}/{f_prob_one_0.denominator} - {f_prob_both_0.numerator}/{f_prob_both_0.denominator} = {f_prob_one_1_one_0.numerator}/{f_prob_one_1_one_0.denominator}\n")

# By symmetry, P(B_x=0, B_y=1) is the same.
f_prob_one_0_one_1 = f_prob_one_1_one_0
print(f"4. By symmetry, P(B_x=0, B_y=1) is also {f_prob_one_0_one_1.numerator}/{f_prob_one_0_one_1.denominator}\n")

# Step 6: Calculate the final probability by summing the results from step 5.
final_prob = f_prob_one_1_one_0 + f_prob_one_0_one_1

print(f"5. The final probability is the sum:")
print(f"   P(B'[i]=1) = {f_prob_one_1_one_0.numerator}/{f_prob_one_1_one_0.denominator} + {f_prob_one_0_one_1.numerator}/{f_prob_one_0_one_1.denominator} = {final_prob.numerator}/{final_prob.denominator}\n")

print(f"The final simplified probability is {final_prob.numerator}/{final_prob.denominator}, which is equal to {float(final_prob)}.")