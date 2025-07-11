from fractions import Fraction

# Step 1: Define the parameters of the Bloom filter
m = 4  # Total bits in the filter B
n = 2  # Number of inputs
k = 2  # Number of hash functions

# Step 2: Calculate the total number of hash operations
N = n * k

print(f"A Bloom filter B has {m} bits, {n} inputs, and {k} hash functions.")
print(f"The total number of hash operations is N = n * k = {N}.")
print("-" * 30)

# Step 3: Calculate the required probabilities using fractions for precision
# Probability that a single hash does NOT land on a specific bit
p_miss_one_bit = Fraction(m - 1, m)
# Probability that a single hash does NOT land on either of two specific bits
p_miss_two_bits = Fraction(m - 2, m)

# Probability that a specific bit remains 0 after all N hashes
p_one_bit_is_0 = p_miss_one_bit ** N
# Probability that two specific bits both remain 0 after all N hashes
p_two_bits_are_0 = p_miss_two_bits ** N

print(f"The analysis for B'[i] = B[i] XOR B[i+2] = 1 requires two bits to be different.")
print(f"The final probability is 2 * (P(one specific bit is 0) - P(two specific bits are 0)).")
print("-" * 30)

# Step 4: Print the intermediate probability calculations
print("Equation for P(one specific bit is 0): ((m-1)/m)^N")
print(f"P(one specific bit is 0) = (({m-1}/{m}))^{N} = {p_one_bit_is_0.numerator}/{p_one_bit_is_0.denominator}")
print()
print("Equation for P(two specific bits are 0): ((m-2)/m)^N")
print(f"P(two specific bits are 0) = (({m-2}/{m}))^{N} = {p_two_bits_are_0.numerator}/{p_two_bits_are_0.denominator}")
print("-" * 30)


# Step 5: Calculate and print the final probability
# P(B'[i]=1) = 2 * (P(B[i]=0) - P(B[i]=0, B[i+2]=0))
final_prob = 2 * (p_one_bit_is_0 - p_two_bits_are_0)

print("Final probability calculation:")
print(f"P(B'[i]=1) = 2 * (P(one bit is 0) - P(two specific bits are 0))")
print(f"P(B'[i]=1) = 2 * ({p_one_bit_is_0.numerator}/{p_one_bit_is_0.denominator} - {p_two_bits_are_0.numerator}/{p_two_bits_are_0.denominator})")
interim_prob = p_one_bit_is_0 - p_two_bits_are_0
print(f"P(B'[i]=1) = 2 * ({interim_prob.numerator}/{interim_prob.denominator})")
print(f"P(B'[i]=1) = {final_prob.numerator}/{final_prob.denominator}")

print("\nThe probability that a randomly chosen bit in B' is 1 is:")
print(final_prob)