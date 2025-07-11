from fractions import Fraction

# The problem reduces to finding three disjoint sets of pairs (i, j)
# for which the condition holds. These sets correspond to three geometric series.

# Set 1: (k, k) for k >= 1.
# The sum of 1/2^(i+j) is the sum of 1/2^(2k) from k=1 to infinity.
# This is a geometric series with a = 1/4 and r = 1/4.
sum1 = Fraction(1, 4) / (1 - Fraction(1, 4))

# Set 2: (k, 2k) for k >= 1.
# The sum of 1/2^(i+j) is the sum of 1/2^(3k) from k=1 to infinity.
# This is a geometric series with a = 1/8 and r = 1/8.
sum2 = Fraction(1, 8) / (1 - Fraction(1, 8))

# Set 3: (2k, k) for k >= 1.
# The sum of 1/2^(i+j) is the sum of 1/2^(3k) from k=1 to infinity.
# This is the same series as for Set 2.
sum3 = Fraction(1, 8) / (1 - Fraction(1, 8))

# The total sum is the sum of the three disjoint series.
total_sum = sum1 + sum2 + sum3

# Output the equation with each number.
s1_num, s1_den = sum1.numerator, sum1.denominator
s2_num, s2_den = sum2.numerator, sum2.denominator
s3_num, s3_den = sum3.numerator, sum3.denominator
total_num, total_den = total_sum.numerator, total_sum.denominator

print(f"The total sum is the sum of three series:")
print(f"Sum 1: {s1_num}/{s1_den}")
print(f"Sum 2: {s2_num}/{s2_den}")
print(f"Sum 3: {s3_num}/{s3_den}")
print(f"Total Sum = {s1_num}/{s1_den} + {s2_num}/{s2_den} + {s3_num}/{s3_den} = {total_num}/{total_den}")
