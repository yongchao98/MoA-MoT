from fractions import Fraction

# Step 1: Define the parameters of the Bloom filter
m = 4  # number of bits in the filter
k = 2  # number of hash functions
n = 2  # number of inputs
total_hashes = k * n

# Step 2: Calculate the probability that a specific bit in B is 0, and then 1
# P(B[i] = 0) = (1 - 1/m)^(k*n)
p_b_is_0 = Fraction(m - 1, m) ** total_hashes

# P(B[i] = 1) = 1 - P(B[i] = 0)
p_b_is_1 = 1 - p_b_is_0

# Step 3: Calculate the joint probability P(B[1]=1, B[3]=1)
# First, find P(B[1]=0 and B[3]=0). This happens if all hashes miss both bits.
# P(B[1]=0, B[3]=0) = (1 - 2/m)^(k*n)
p_b1_and_b3_are_0 = Fraction(m - 2, m) ** total_hashes

# Next, find P(B[1]=0 or B[3]=0) using inclusion-exclusion.
# P(A or C) = P(A) + P(C) - P(A and C)
p_b1_or_b3_are_0 = p_b_is_0 + p_b_is_0 - p_b1_and_b3_are_0

# Finally, P(B[1]=1, B[3]=1) is the complement of P(B[1]=0 or B[3]=0).
p_b1_and_b3_are_1 = 1 - p_b1_or_b3_are_0

# Step 4: Calculate the final probability P(B'[i] = 1)
# P(B'[i] = 1) = P(B[1] != B[3]) = P(B[1]=1) + P(B[3]=1) - 2 * P(B[1]=1, B[3]=1)
# By symmetry, P(B[1]=1) = P(B[3]=1) = p_b_is_1
final_prob = p_b_is_1 + p_b_is_1 - 2 * p_b1_and_b3_are_1

# Print the final equation with all intermediate values
# To avoid clutter, intermediate results are shown in their simplified fractional form
equation = (
    f"P(B'[i]=1) = P(B[1]=1) + P(B[3]=1) - 2 * P(B[1]=1, B[3]=1)\n"
    f"           = {p_b_is_1} + {p_b_is_1} - 2 * {p_b1_and_b3_are_1}\n"
    f"           = {p_b_is_1 + p_b_is_1} - {2 * p_b1_and_b3_are_1}\n"
    f"           = {final_prob}"
)

print("The calculation is as follows:")
print(equation)
print("\nThe final probability is:")
print(float(final_prob))
print(f"Which is equivalent to the fraction: {final_prob}")
