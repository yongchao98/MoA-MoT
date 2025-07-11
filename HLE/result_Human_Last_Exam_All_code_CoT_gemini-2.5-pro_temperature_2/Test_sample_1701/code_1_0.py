import math

# --- Parameters of the Bloom Filter ---
m = 4  # Total bits in the Bloom filter B
k = 2  # Number of hash functions
n = 2  # Number of inputs

# --- Step 1: Calculate the probability of a bit in B being 0 or 1 ---

# The probability that a bit remains 0 after all insertions is calculated
# based on the probability that no hash function maps to its position.
# We will calculate p0 = P(B[j] = 0) as an exact fraction to avoid rounding errors.
p0_num = (m - 1)**(k * n)
p0_den = m**(k * n)

# The probability that a bit is 1 is p1 = P(B[j] = 1) = 1 - p0.
p1_num = p0_den - p0_num
p1_den = p0_den

# --- Step 2: Calculate the probability of a bit in B' being 1 ---

# A bit B'[i] is 1 if B_1[i] and B_2[i] are different.
# P(B'[i]=1) = P(B_1[i]=0, B_2[i]=1) + P(B_1[i]=1, B_2[i]=0)
# Assuming independence of the bits in B, this simplifies to:
# P(B'[i]=1) = 2 * P(B[j]=0) * P(B[j]=1)

# --- Step 3: Compute the final numerical answer ---

# Calculate the final probability using the fractions derived above.
final_prob_num = 2 * p0_num * p1_num
final_prob_den = p0_den * p1_den

# For clarity, simplify the resulting fraction by finding the greatest common divisor.
common_divisor = math.gcd(final_prob_num, final_prob_den)
simplified_num = final_prob_num // common_divisor
simplified_den = final_prob_den // common_divisor

# Also calculate the final decimal value.
final_prob_decimal = final_prob_num / final_prob_den

# --- Step 4: Print the full calculation and result ---

print("Calculating the probability that a random bit B'[i] is 1.")
print("-" * 50)
print(f"Given parameters: m={m}, k={k}, n={n}\n")

print("First, we find the probability of a single bit in the original Bloom filter B being 0 (p0) or 1 (p1).")
print(f"p0 = P(B[j]=0) = (1 - 1/{m})^({k}*{n}) = ({m - 1}/{m})^{k*n} = {p0_num}/{p0_den}")
print(f"p1 = P(B[j]=1) = 1 - p0 = {p1_num}/{p1_den}\n")

print("Next, we find the probability of a bit in B' being 1.")
print("P(B'[i]=1) = P(B_1[i] != B_2[i]) = 2 * p0 * p1\n")

print("Final Equation with the numbers inserted:")
print(f"P(B'[i]=1) = 2 * ({p0_num}/{p0_den}) * ({p1_num}/{p1_den})")
print(f"           = {final_prob_num}/{final_prob_den}")
print(f"           = {simplified_num}/{simplified_den} (in simplest form)")
print(f"           = {final_prob_decimal} (as a decimal)\n")

print(f"The final probability is {final_prob_decimal:.5f}.")