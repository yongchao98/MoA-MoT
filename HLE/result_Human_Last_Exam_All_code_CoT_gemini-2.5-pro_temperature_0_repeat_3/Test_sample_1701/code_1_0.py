import fractions

# Problem parameters
m = 4  # Number of bits in the Bloom filter B
k = 2  # Number of hash functions
n = 2  # Number of inputs

# Total number of hash operations
total_hashes = n * k

# The goal is to find P(B'[i] = 1). Due to symmetry, we can calculate P(B'[1] = 1).
# B'[1] = B[1] XOR B[3]. This is 1 if and only if B[1] != B[3].
# P(B'[1] = 1) = P(B[1] != B[3]) = P(B[1]=0, B[3]=1) + P(B[1]=1, B[3]=0)

# Step 1: Calculate the probability that both B[1] and B[3] are 0.
# For a single hash, the probability of not hitting index 1 or 3 is (m-2)/m.
# For all hashes, this is ((m-2)/m) ^ total_hashes.
p_b1_0_and_b3_0_frac = fractions.Fraction(m - 2, m) ** total_hashes

# Step 2: Calculate the probability that B[1] is 0.
# For a single hash, the probability of not hitting index 1 is (m-1)/m.
# For all hashes, this is ((m-1)/m) ^ total_hashes.
p_b1_0_frac = fractions.Fraction(m - 1, m) ** total_hashes

# Step 3: Calculate P(B[1]=0, B[3]=1)
# We know P(B[1]=0) = P(B[1]=0, B[3]=0) + P(B[1]=0, B[3]=1).
# So, P(B[1]=0, B[3]=1) = P(B[1]=0) - P(B[1]=0, B[3]=0).
p_b1_0_and_b3_1_frac = p_b1_0_frac - p_b1_0_and_b3_0_frac

# Step 4: By symmetry, P(B[1]=1, B[3]=0) is the same.
p_b1_1_and_b3_0_frac = p_b1_0_and_b3_1_frac

# Step 5: Calculate the final probability P(B[1] != B[3]).
final_prob_frac = p_b1_0_and_b3_1_frac + p_b1_1_and_b3_0_frac

# Print the derivation with the calculated numbers
print("The probability of a randomly chosen bit in B' being 1 is P(B'[i]=1).")
print("Let's calculate this for i=1: P(B'[1]=1) = P(B[1] != B[3])")
print("P(B[1] != B[3]) = P(B[1]=0, B[3]=1) + P(B[1]=1, B[3]=0)")
print("\nFirst, we find the probability of a bit (or bits) being 0 after all insertions.")
print(f"P(B[1]=0) = ((m-1)/m)^(n*k) = (({m}-1)/{m})^({n}*{k}) = {p_b1_0_frac.numerator}/{p_b1_0_frac.denominator}")
print(f"P(B[1]=0, B[3]=0) = ((m-2)/m)^(n*k) = (({m}-2)/{m})^({n}*{k}) = {p_b1_0_and_b3_0_frac.numerator}/{p_b1_0_and_b3_0_frac.denominator}")
print("\nNow we can calculate the required joint probabilities:")
print(f"P(B[1]=0, B[3]=1) = P(B[1]=0) - P(B[1]=0, B[3]=0) = {p_b1_0_frac} - {p_b1_0_and_b3_0_frac} = {p_b1_0_and_b3_1_frac}")
print(f"By symmetry, P(B[1]=1, B[3]=0) = {p_b1_1_and_b3_0_frac}")
print("\nFinally, we sum these probabilities for the final answer:")
print(f"P(B'[1]=1) = {p_b1_0_and_b3_1_frac} + {p_b1_1_and_b3_0_frac} = {final_prob_frac}")
print(f"\nThe final probability is {final_prob_frac.numerator}/{final_prob_frac.denominator}, which is approximately {float(final_prob_frac):.5f}.")
