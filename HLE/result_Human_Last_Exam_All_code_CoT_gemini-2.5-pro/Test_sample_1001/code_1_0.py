from fractions import Fraction

# The problem reduces to finding three families of pairs (i, j) that satisfy the condition.
# The sum is calculated over these three disjoint families.

# Family 1: (d, d) for d >= 1.
# The sum is sum_{d=1 to inf} 1 / 2^(d+d) = sum_{d=1 to inf} (1/4)^d.
# This is a geometric series with a=1/4, r=1/4.
a1 = Fraction(1, 4)
r1 = Fraction(1, 4)
sum1 = a1 / (1 - r1)

# Family 2: (d, 2d) for d >= 1.
# The sum is sum_{d=1 to inf} 1 / 2^(d+2d) = sum_{d=1 to inf} (1/8)^d.
# This is a geometric series with a=1/8, r=1/8.
a2 = Fraction(1, 8)
r2 = Fraction(1, 8)
sum2 = a2 / (1 - r2)

# Family 3: (2d, d) for d >= 1.
# The sum is the same as for Family 2.
sum3 = sum2

# The total sum is the sum of the three series.
total_sum = sum1 + sum2 + sum3

# As requested, we print each number in the final equation.
print("The total sum is the sum of three geometric series corresponding to the three families of pairs.")
print(f"Sum = {sum1.numerator}/{sum1.denominator} + {sum2.numerator}/{sum2.denominator} + {sum3.numerator}/{sum3.denominator}")
print(f"Result = {total_sum.numerator}/{total_sum.denominator}")
