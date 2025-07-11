import math

# Define the parameters of the Bloom filter
m = 4  # Total bits in B
k = 2  # Number of hash functions
n = 2  # Number of inputs

# Total number of hash operations
nk = n * k

# Calculate P(B[j]=0), the probability that a single specific bit in B is 0.
# This happens if none of the nk hashes land on this bit.
# The probability for one hash to not land on the bit is (m-1)/m.
p_bit_is_0 = ((m - 1) / m) ** nk
p_bit_is_0_num = (m - 1) ** nk
p_bit_is_0_den = m ** nk


# Calculate P(B[j]=0, B[l]=0), the probability that two distinct specific bits in B are both 0.
# This happens if none of the nk hashes land on either of these two bits.
# The probability for one hash not to land on either bit is (m-2)/m.
p_two_bits_are_0 = ((m - 2) / m) ** nk
p_two_bits_are_0_num = (m - 2) ** nk
p_two_bits_are_0_den = m ** nk


# The probability that B'[i] = 1 is P(B[j] != B[l]).
# This is P(B[j]=1, B[l]=0) + P(B[j]=0, B[l]=1).
# By symmetry, this is 2 * P(B[j]=1, B[l]=0).
# We also know P(B[l]=0) = P(B[j]=1, B[l]=0) + P(B[j]=0, B[l]=0).
# So, P(B[j]=1, B[l]=0) = P(B[l]=0) - P(B[j]=0, B[l]=0).
# The final probability is 2 * (p_bit_is_0 - p_two_bits_are_0).
final_prob = 2 * (p_bit_is_0 - p_two_bits_are_0)
final_prob_num = 2 * (p_bit_is_0_num * (p_two_bits_are_0_den // p_bit_is_0_den) - p_two_bits_are_0_num)
final_prob_den = p_two_bits_are_0_den
common_divisor = math.gcd(final_prob_num, final_prob_den)
final_prob_frac_num = final_prob_num // common_divisor
final_prob_frac_den = final_prob_den // common_divisor


print(f"Given parameters: m={m}, k={k}, n={n}")
print(f"Total hashes (n*k) = {nk}\n")
print("Let's calculate P(B'[i]=1) = P(one bit is 1) - P(two bits are 0)")
print("The probability is given by the formula:")
print(f"P(B'[i]=1) = 2 * (P(B[j]=0) - P(B[j]=0, B[l]=0))")
print(f"P(B'[i]=1) = 2 * ( (({m}-1)/{m})^({n}*{k}) - (({m}-2)/{m})^({n}*{k}) )")
print(f"P(B'[i]=1) = 2 * ( ({m-1}/{m})^{nk} - ({m-2}/{m})^{nk} )")
print(f"P(B'[i]=1) = 2 * ( (({p_bit_is_0_num}/{p_bit_is_0_den})) - (({p_two_bits_are_0_num}/{p_two_bits_are_0_den})) )")
print(f"P(B'[i]=1) = 2 * ( {p_bit_is_0_num * (p_two_bits_are_0_den // p_bit_is_0_den) - p_two_bits_are_0_num}/{p_two_bits_are_0_den} )")
print(f"P(B'[i]=1) = {final_prob_num}/{final_prob_den} = {final_prob_frac_num}/{final_prob_frac_den}\n")
print(f"The final probability is: {final_prob}")
