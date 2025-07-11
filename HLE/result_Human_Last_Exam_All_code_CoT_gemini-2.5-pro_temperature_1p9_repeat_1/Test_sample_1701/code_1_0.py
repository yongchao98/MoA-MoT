import math

# Step 1: Define System Parameters
m = 4  # Total bits in the Bloom filter B
k = 2  # Number of hash functions
n = 2  # Number of inputs

# Total number of hash computations
N = n * k

# Step 2 & 3: We want to find P(B'[i]=1) = 2 * P(B_1[i]=1, B_2[i]=0)
# Step 4: This is equal to 2 * (P(B_2[i]=0) - P(B_1[i]=0, B_2[i]=0))

# Step 5: Calculate the fundamental probabilities using their numerators and denominators

# Probability that a single bit remains 0: P(B_2[i]=0) = ((m-1)/m)^N
p_zero_num = (m - 1)**N
p_zero_den = m**N

# Probability that two specific bits remain 0: P(B_1[i]=0, B_2[i]=0) = ((m-2)/m)^N
p_zero_zero_num = (m - 2)**N
p_zero_zero_den = m**N

# Step 6: Combine and solve
# P(B_1[i]=1, B_2[i]=0) = P(B_2[i]=0) - P(B_1[i]=0, B_2[i]=0)
p_one_zero_num = p_zero_num - p_zero_zero_num
p_one_zero_den = p_zero_den

# The final probability is 2 * P(B_1[i]=1, B_2[i]=0)
final_prob_num = 2 * p_one_zero_num
final_prob_den = p_one_zero_den

# Simplify the final fraction
common_divisor = math.gcd(final_prob_num, final_prob_den)
s_final_prob_num = final_prob_num // common_divisor
s_final_prob_den = final_prob_den // common_divisor

# Print the calculation steps and the final equation
print("The probability P(B'[i]=1) is calculated as 2 * (P(a specific bit = 0) - P(two specific bits = 0)).")
print("\nIntermediate probabilities:")
print(f"P(a specific bit = 0) = (({m}-1)/{m})^({n}*{k}) = ({m - 1}^{N}) / ({m}^{N}) = {p_zero_num}/{p_zero_den}")
print(f"P(two specific bits = 0) = (({m}-2)/{m})^({n}*{k}) = ({m - 2}^{N}) / ({m}^{N}) = {p_zero_zero_num}/{p_zero_zero_den}")
print("\nFinal calculation:")
print(f"P(B'[i]=1) = 2 * ({p_zero_num}/{p_zero_den} - {p_zero_zero_num}/{p_zero_zero_den}) = 2 * ({p_one_zero_num}/{p_one_zero_den}) = {final_prob_num}/{final_prob_den}")
print(f"\nThe simplified probability is: {s_final_prob_num}/{s_final_prob_den}")
