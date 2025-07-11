from fractions import Fraction

# The problem asks for the sum of 1/2^(i+j) over a set S of pairs (i, j).
# The derivation shows that the set S consists of three disjoint sets of pairs:
# 1. Pairs of the form (k, k) for k >= 1.
# 2. Pairs of the form (k, 2k) for k >= 1.
# 3. Pairs of the form (2k, k) for k >= 1.
# We calculate the sum for each of these sets and add them together.

# Part 1: Sum for pairs (k, k)
# The sum is Sum_{k=1 to inf} 1/2^(k+k) = Sum_{k=1 to inf} (1/4)^k.
# This is a geometric series with first term a = 1/4 and common ratio r = 1/4.
# The sum is a / (1 - r).
sum_part1 = Fraction(1, 4) / (1 - Fraction(1, 4))

# Part 2: Sum for pairs (k, 2k)
# The sum is Sum_{k=1 to inf} 1/2^(k+2k) = Sum_{k=1 to inf} (1/8)^k.
# This is a geometric series with first term a = 1/8 and common ratio r = 1/8.
# The sum is a / (1 - r).
sum_part2 = Fraction(1, 8) / (1 - Fraction(1, 8))

# Part 3: Sum for pairs (2k, k)
# The sum is Sum_{k=1 to inf} 1/2^(2k+k) = Sum_{k=1 to inf} (1/8)^k.
# This is identical to the sum for Part 2.
sum_part3 = sum_part2

# The total sum is the sum of these three parts, as the sets of pairs are disjoint.
total_sum = sum_part1 + sum_part2 + sum_part3

# Final output: an equation showing the sum of the components.
print(f"The total sum is the sum of three parts: {sum_part1} + {sum_part2} + {sum_part3} = {total_sum}")
