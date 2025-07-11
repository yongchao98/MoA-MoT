import math

# Problem parameters
m = 4  # number of bits in B
k = 2  # number of hash functions
n = 2  # number of inputs
N = k * n # total number of hashings

print(f"Step 1: The goal is to find P(B'[i]=1). We'll compute P(B[1] xor B[3] = 1).")
print(f"This is equal to P(B[1]=1, B[3]=0) + P(B[1]=0, B[3]=1).\n")

print(f"Step 2: Calculate P(B[1]=1, B[3]=0).")
print(f"This is P(no hash hits bit 3) - P(no hash hits bit 1 or 3).\n")

# P(no hash hits bit 3) = (1 - 1/m)^N
prob_no_hit_3_num = (m - 1)**N
prob_no_hit_3_den = m**N
print(f"P(no hash hits bit 3) = (1 - 1/{m})^{N} = ({(m-1)}/{m})^{N} = {prob_no_hit_3_num}/{prob_no_hit_3_den}")

# P(no hash hits bit 1 or 3) = (1 - 2/m)^N
# We simplify the base fraction (m-2)/m = 2/4 = 1/2 first
base_num = m-2
base_den = m
common = math.gcd(base_num, base_den)
prob_no_hit_1_or_3_base_num = base_num // common
prob_no_hit_1_or_3_base_den = base_den // common

prob_no_hit_1_or_3_num = prob_no_hit_1_or_3_base_num**N
prob_no_hit_1_or_3_den = prob_no_hit_1_or_3_base_den**N
print(f"P(no hash hits bit 1 or 3) = (1 - 2/{m})^{N} = ({prob_no_hit_1_or_3_base_num}/{prob_no_hit_1_or_3_base_den})^{N} = {prob_no_hit_1_or_3_num}/{prob_no_hit_1_or_3_den}\n")

# P(b1=1, b3=0) = prob_no_hit_3 - prob_no_hit_1_or_3
# To subtract, we find a common denominator
common_den = prob_no_hit_3_den # This is lcm(256, 16) = 256
term2_num_scaled = prob_no_hit_1_or_3_num * (common_den // prob_no_hit_1_or_3_den)
joint_prob_num = prob_no_hit_3_num - term2_num_scaled
joint_prob_den = common_den
print(f"P(B[1]=1, B[3]=0) = {prob_no_hit_3_num}/{prob_no_hit_3_den} - {prob_no_hit_1_or_3_num}/{prob_no_hit_1_or_3_den} = {prob_no_hit_3_num}/{common_den} - {term2_num_scaled}/{common_den} = {joint_prob_num}/{joint_prob_den}\n")

# Final Probability = 2 * P(b1=1, b3=0)
final_prob_num = 2 * joint_prob_num
final_prob_den = joint_prob_den

# Simplify final fraction
final_common_divisor = math.gcd(final_prob_num, final_prob_den)
final_simp_num = final_prob_num // final_common_divisor
final_simp_den = final_prob_den // final_common_divisor

print(f"Step 3: Sum the probabilities.")
print(f"P(B'[1]=1) = P(B[1]=1, B[3]=0) + P(B[1]=0, B[3]=1)")
print(f"= 2 * P(B[1]=1, B[3]=0)")
print(f"= 2 * {joint_prob_num}/{joint_prob_den}")
print(f"= {final_prob_num}/{final_prob_den}\n")

print(f"Final Answer:")
print(f"The final simplified equation for the probability is:")
print(f"P(B'[i]=1) = 2 * [ (1 - 1/{m})^{k*n} - (1 - 2/{m})^{k*n} ] = 2 * [ ({(m-1)}/{m})^{N} - ({(m-2)}/{m})^{N} ] = {final_simp_num}/{final_simp_den}")
print(f"The decimal probability is: {final_prob_num/final_prob_den}")
