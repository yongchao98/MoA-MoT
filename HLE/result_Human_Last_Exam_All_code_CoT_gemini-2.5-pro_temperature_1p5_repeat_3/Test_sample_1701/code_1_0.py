import fractions

# Parameters from the problem description
m = 4  # number of bits in the Bloom filter B
k = 2  # number of hash functions
n = 2  # number of inputs

# --- Step 1: Calculate the probability of a bit being 0 (p_0) ---
# The probability that a specific bit is 0 is (1 - 1/m)^(k*n).
# We use the fractions module for exact rational number arithmetic.
p_0 = fractions.Fraction(m - 1, m) ** (k * n)

# --- Step 2: Calculate the probability of a bit being 1 (p_1) ---
p_1 = 1 - p_0

# --- Step 3: Calculate the probability of a bit in B' being 1 ---
# B'[i] = B1[i] XOR B2[i]. This is 1 if and only if B1[i] != B2[i].
# The probability is P(B1[i]=0, B2[i]=1) + P(B1[i]=1, B2[i]=0).
# Assuming independence of the bits in B, this is p_0*p_1 + p_1*p_0 = 2*p_0*p_1.
final_prob = 2 * p_0 * p_1

# --- Step 4: Output the results as requested ---
print(f"The probability of a bit in the Bloom filter B being 0 is p_0.")
print(f"p_0 = (1 - 1/{m})^({k}*{n}) = {p_0.numerator}/{p_0.denominator}")
print("-" * 20)
print(f"The probability of a bit in the Bloom filter B being 1 is p_1.")
print(f"p_1 = 1 - p_0 = {p_1.numerator}/{p_1.denominator}")
print("-" * 20)
print("The probability of a bit in B' being 1 is P(B'[i]=1) = 2 * p_0 * p_1.")
print("The final equation with the calculated numbers is:")
print(f"P(B'[i]=1) = 2 * ({p_0.numerator}/{p_0.denominator}) * ({p_1.numerator}/{p_1.denominator}) = {final_prob.numerator}/{final_prob.denominator}")
print("-" * 20)
print(f"The final probability is {final_prob.numerator}/{final_prob.denominator}")
print(f"As a decimal, this is: {float(final_prob)}")
