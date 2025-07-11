from fractions import Fraction

# The set S consists of three types of pairs (i, j):
# 1. (k, k) for k >= 1
# 2. (k, 2k) for k >= 1
# 3. (2k, k) for k >= 1
# We calculate the sum of 1/2^(i+j) for each set and add them up.

# Part 1: Sum for pairs (k, k)
# This corresponds to the sum of 1/2^(2k) from k=1 to infinity.
# It is a geometric series with a = 1/4 and r = 1/4.
a1 = Fraction(1, 4)
r1 = Fraction(1, 4)
sum1 = a1 / (1 - r1)

# Part 2: Sum for pairs (k, 2k)
# This corresponds to the sum of 1/2^(3k) from k=1 to infinity.
# It is a geometric series with a = 1/8 and r = 1/8.
a2 = Fraction(1, 8)
r2 = Fraction(1, 8)
sum2 = a2 / (1 - r2)

# Part 3: Sum for pairs (2k, k)
# This corresponds to the sum of 1/2^(3k) from k=1 to infinity.
# It is a geometric series with a = 1/8 and r = 1/8.
a3 = Fraction(1, 8)
r3 = Fraction(1, 8)
sum3 = a3 / (1 - r3)

# The total sum is the sum of the three disjoint parts.
total_sum = sum1 + sum2 + sum3

print("The calculation is the sum of three geometric series:")
print(f"Sum for pairs (k,k): {sum1}")
print(f"Sum for pairs (k,2k): {sum2}")
print(f"Sum for pairs (2k,k): {sum3}")
print("The final equation with each number is:")
print(f"{sum1} + {sum2} + {sum3} = {total_sum}")
