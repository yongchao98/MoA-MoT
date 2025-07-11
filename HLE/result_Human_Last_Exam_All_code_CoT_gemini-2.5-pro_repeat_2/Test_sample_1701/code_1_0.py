from fractions import Fraction

# Problem parameters
m = 4  # bits in the Bloom filter
k = 2  # hash functions
n = 2  # inputs

# The total number of hash computations
kn = k * n

# The probability that a specific bit in B is 0.
# This happens if none of the kn hashes point to this bit.
# The probability for a single hash NOT to point to this bit is (m-1)/m.
# For kn hashes, it's ((m-1)/m)**kn.
p0_num = (m - 1)**kn
p0_den = m**kn
p0 = Fraction(p0_num, p0_den)

# The probability that two specific bits in B are both 0.
# This happens if none of the kn hashes point to either of these two bits.
# The probability for a single hash NOT to point to either bit is (m-2)/m.
# For kn hashes, it's ((m-2)/m)**kn.
p00_num = (m - 2)**kn
p00_den = m**kn
p00 = Fraction(p00_num, p00_den)

# We want to find P(B'[i]=1). Let's analyze for i=1.
# P(B'[1]=1) = P(B[1] XOR B[3] = 1)
# This is the probability that B[1] and B[3] are different.
# P(B[1] != B[3]) = P(B[1]=0, B[3]=1) + P(B[1]=1, B[3]=0)

# Let's find P(B[1]=0, B[3]=1).
# We know P(B[1]=0) = P(B[1]=0, B[3]=0) + P(B[1]=0, B[3]=1).
# So, P(B[1]=0, B[3]=1) = P(B[1]=0) - P(B[1]=0, B[3]=0)
p01 = p0 - p00

# By symmetry, P(B[1]=1, B[3]=0) is the same as P(B[1]=0, B[3]=1).
p10 = p01

# The final probability is the sum.
final_prob = p01 + p10

# Print the explanation and the final equation
print("The goal is to find the probability P(B'[i] = 1).")
print("By construction, B'[i] = B_1[i] XOR B_2[i].")
print("For i=1, this means B'[1] = B[1] XOR B[3].")
print("P(B'[1] = 1) is the probability that B[1] and B[3] have different values.")
print("This can be expressed as P(B[1] != B[3]) = P(B[1]=0, B[3]=1) + P(B[1]=1, B[3]=0).\n")

print(f"Let p0 = P(B[j]=0) and p00 = P(B[j]=0, B[l]=0) for j!=l.")
print(f"With m={m}, k={k}, n={n}:")
print(f"p0 = (1 - 1/{m})^({k}*{n}) = ({m-1}/{m})^{kn} = {p0_num}/{p0_den}")
print(f"p00 = (1 - 2/{m})^({k}*{n}) = ({m-2}/{m})^{kn} = {p00_num}/{p00_den}\n")

print(f"From these, we can find the probability of one bit being 0 and the other 1:")
print(f"P(B[1]=0, B[3]=1) = p0 - p00 = {p0} - {p00} = {p01}\n")

print(f"By symmetry, P(B[1]=1, B[3]=0) is the same.")
print(f"The final probability is P(B[1] != B[3]) = 2 * P(B[1]=0, B[3]=1)")
print(f"P(B'[i]=1) = 2 * {p01.numerator}/{p01.denominator} = {final_prob.numerator}/{final_prob.denominator}")

# Final numerical answer
print("\nThe final probability as a decimal is:")
print(float(final_prob))