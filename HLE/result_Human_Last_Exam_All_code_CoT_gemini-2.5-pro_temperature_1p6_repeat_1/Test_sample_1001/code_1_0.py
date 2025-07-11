import fractions

# Based on the analysis, the set S consists of three disjoint types of pairs.
# We calculate the sum of 1/2^(i+j) for each type and add them up.

# Type 1: Pairs (k, k) for k >= 1.
# The sum is a geometric series: sum_{k=1 to inf} (1/4)^k.
# The sum of a geometric series a + ar + ar^2 + ... is a / (1 - r).
# Here, a = 1/4 and r = 1/4.
sum1 = fractions.Fraction(1, 4) / (1 - fractions.Fraction(1, 4))

# Type 2: Pairs (k, 2k) for k >= 1.
# The sum is a geometric series: sum_{k=1 to inf} (1/8)^k.
# Here, a = 1/8 and r = 1/8.
sum2 = fractions.Fraction(1, 8) / (1 - fractions.Fraction(1, 8))

# Type 3: Pairs (2k, k) for k >= 1.
# The sum is the same as for Type 2.
sum3 = sum2

# The total sum is the sum of the sums for the three disjoint sets.
total_sum = sum1 + sum2 + sum3

# We output the numbers that form the final equation.
# The three sums are 1/3, 1/7, and 1/7.
# Their total sum is 1/3 + 1/7 + 1/7 = 7/21 + 6/21 = 13/21.
print(f"{sum1.numerator}/{sum1.denominator} + {sum2.numerator}/{sum2.denominator} + {sum3.numerator}/{sum3.denominator} = {total_sum.numerator}/{total_sum.denominator}")