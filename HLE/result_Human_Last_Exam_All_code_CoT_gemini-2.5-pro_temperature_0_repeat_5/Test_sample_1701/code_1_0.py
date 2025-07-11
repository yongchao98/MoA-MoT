from fractions import Fraction

# Define the parameters of the Bloom filter as given in the problem.
m = 4  # Total bits in the original filter B
k = 2  # Number of hash functions
n = 2  # Number of inputs

# Step 1: Calculate the probability that a specific bit in the original Bloom filter B is 0.
# A bit is 0 only if none of the n*k hashings point to it.
# The probability of a single hash NOT pointing to a specific bit is (1 - 1/m).
# With n*k total hashings, this probability is (1 - 1/m)^(n*k).
# We use the Fraction class for precise arithmetic.
p0 = Fraction(1 - 1/m)**(n * k)

# Step 2: Calculate the probability that a specific bit in B is 1.
# This is the complement of the probability of it being 0.
p1 = 1 - p0

# Step 3: Calculate the probability that a bit in the derived filter B' is 1.
# A bit B'[i] is 1 if the two corresponding bits from B1 and B2 are different.
# P(B'[i] = 1) = P(bit_B1=0, bit_B2=1) + P(bit_B1=1, bit_B2=0)
# Assuming independence of bits in the original filter, this becomes:
# P(B'[i] = 1) = p0 * p1 + p1 * p0 = 2 * p0 * p1
final_prob = 2 * p0 * p1

# Step 4: Print the final equation with all the numbers, as requested.
# We show the fractional values for clarity and precision.
print(f"The probability of a bit in the original filter B being 0 is p0 = (1 - 1/{m})^({n}*{k}) = {p0.numerator}/{p0.denominator}")
print(f"The probability of a bit in the original filter B being 1 is p1 = 1 - p0 = {p1.numerator}/{p1.denominator}")
print("\nThe probability of a bit in the derived filter B' being 1 is given by the equation: 2 * p0 * p1")
print("\nSubstituting the fractional values, the final equation is:")
print(f"2 * ({p0.numerator}/{p0.denominator}) * ({p1.numerator}/{p1.denominator}) = {final_prob.numerator}/{final_prob.denominator}")

# Also print the decimal representation for context.
print(f"\nThe final probability as a decimal is approximately: {float(final_prob)}")