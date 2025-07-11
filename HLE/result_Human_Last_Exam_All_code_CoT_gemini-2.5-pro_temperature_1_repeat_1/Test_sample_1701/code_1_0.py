from fractions import Fraction

# Step 1: Define the parameters for the Bloom filter.
m = 4  # Number of bits in the filter
k = 2  # Number of hash functions
n = 2  # Number of inputs inserted

# Total number of hash computations
nk = n * k

# Step 2: Calculate the probability that a specific bit B[j] is 0.
# A bit is 0 if none of the nk hashes land on it.
# The probability for one hash to miss a specific bit is (1 - 1/m).
# For nk hashes, this is P(B[j]=0) = (1 - 1/m)^(n*k).
p0 = Fraction((m - 1)**nk, m**nk)

# Step 3: Calculate the probability that two distinct bits, B[i] and B[j], are both 0.
# A hash misses both bits with probability (1 - 2/m).
# For nk hashes, this is P(B[i]=0, B[j]=0) = (1 - 2/m)^(n*k).
p00 = Fraction((m - 2)**nk, m**nk)

# Step 4: Calculate the probability for the derived bit B'[i] to be 1.
# B'[i] = 1 means B[i] XOR B[j] = 1, which happens if (B[i]=0, B[j]=1) or (B[i]=1, B[j]=0).
# From the law of total probability, P(B[i]=0) = P(B[i]=0, B[j]=0) + P(B[i]=0, B[j]=1).
# Therefore, P(B[i]=0, B[j]=1) = P(B[i]=0) - P(B[i]=0, B[j]=0) = p0 - p00.
# By symmetry, P(B[i]=1, B[j]=0) is the same.
# The final probability is P(B'[i]=1) = 2 * (p0 - p00).
final_prob = 2 * (p0 - p00)

# Step 5: Print the final equation with all the numbers.
print(f"The probability that a randomly chosen bit in B' is 1 can be calculated as follows:")
print(f"P(B'[i]=1) = 2 * (P(B[j]=0) - P(B[i]=0, B[j]=0))")
print(f"P(B[j]=0) = (1 - 1/{m})^({n*k}) = {p0.numerator}/{p0.denominator}")
print(f"P(B[i]=0, B[j]=0) = (1 - 2/{m})^({n*k}) = {p00.numerator}/{p00.denominator}")
print("\nFinal Equation:")
print(f"P(B'[i]=1) = 2 * ({p0.numerator}/{p0.denominator} - {p00.numerator}/{p00.denominator}) = {final_prob.numerator}/{final_prob.denominator}")

# Print the final result in the requested format
print(f"\nThe calculated probability is {final_prob}.")
print(f"\n<<<65/128>>>")