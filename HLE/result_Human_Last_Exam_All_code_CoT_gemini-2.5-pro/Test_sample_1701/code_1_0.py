import fractions

# Parameters of the Bloom Filter
m = 4  # number of bits
k = 2  # number of hash functions
n = 2  # number of inputs

# Total number of hash computations
nk = n * k

# --- Step 1: Calculate P(B[i]=0) ---
# Probability that a single hash does not land on a specific bit
prob_single_hash_miss = (m - 1) / m
# Probability that all nk hashes do not land on a specific bit
prob_b_i_is_0 = prob_single_hash_miss ** nk

# --- Step 2: Calculate P(B[i]=0, B[j]=0) for i != j ---
# Probability that a single hash does not land on two specific bits
prob_single_hash_miss_two = (m - 2) / m
# Probability that all nk hashes do not land on two specific bits
prob_b_i_and_b_j_are_0 = prob_single_hash_miss_two ** nk

# --- Step 3: Calculate P(B[i]=0, B[j]=1) ---
# P(B[i]=0, B[j]=1) = P(B[i]=0) - P(B[i]=0, B[j]=0)
prob_b_i_is_0_and_b_j_is_1 = prob_b_i_is_0 - prob_b_i_and_b_j_are_0

# --- Step 4: Calculate the final probability P(B'[i]=1) ---
# P(B'[i]=1) = P(B[i] xor B[j] = 1) = 2 * P(B[i]=0, B[j]=1)
final_prob = 2 * prob_b_i_is_0_and_b_j_is_1

# --- Print the explanation and results ---
print("The problem is to find P(B'[i]=1), which is P(B[i] xor B[j] = 1) for corresponding bits i and j.")
print("The formula is: P(B'[i]=1) = 2 * (P(B[i]=0) - P(B[i]=0, B[j]=0))\n")

print("1. Calculate P(B[i]=0):")
print(f"   P(B[i]=0) = (({m}-1)/{m})^({n}*{k}) = ({fractions.Fraction(m - 1, m)})^4 = {fractions.Fraction(prob_b_i_is_0)} = {prob_b_i_is_0}")

print("\n2. Calculate P(B[i]=0, B[j]=0):")
print(f"   P(B[i]=0, B[j]=0) = (({m}-2)/{m})^({n}*{k}) = ({fractions.Fraction(m - 2, m)})^4 = {fractions.Fraction(prob_b_i_and_b_j_are_0)} = {prob_b_i_and_b_j_are_0}")

print("\n3. Calculate the final probability:")
print(f"   P(B'[i]=1) = 2 * (P(B[i]=0) - P(B[i]=0, B[j]=0))")
print(f"   P(B'[i]=1) = 2 * ({prob_b_i_is_0} - {prob_b_i_and_b_j_are_0})")
print(f"   P(B'[i]=1) = 2 * ({prob_b_i_is_0_and_b_j_is_1})")
print(f"   P(B'[i]=1) = {final_prob}")
print(f"\nThe final probability as a fraction is: {fractions.Fraction(final_prob).limit_denominator()}")
