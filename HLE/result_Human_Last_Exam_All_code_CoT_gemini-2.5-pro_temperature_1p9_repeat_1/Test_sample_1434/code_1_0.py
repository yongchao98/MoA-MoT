from fractions import Fraction

# The problem boils down to calculating the sum of three geometric series.
# The sum is given by: S = 1/2 * Sum_1 + 1/2 * Sum_2 + Sum_3

# Sum_1 is the sum of (3/4)^n from n=1 to infinity.
# The sum of a geometric series a + ar + ar^2 + ... is a/(1-r).
# Here, a = 3/4 and r = 3/4.
a1 = Fraction(3, 4)
r1 = Fraction(3, 4)
sum1 = a1 / (1 - r1)

# Sum_2 is the sum of (-1/4)^n from n=1 to infinity.
# Here, a = -1/4 and r = -1/4.
a2 = Fraction(-1, 4)
r2 = Fraction(-1, 4)
sum2 = a2 / (1 - r2)

# Sum_3 is the sum of (1/4)^n from n=1 to infinity.
# Here, a = 1/4 and r = 1/4.
a3 = Fraction(1, 4)
r3 = Fraction(1, 4)
sum3 = a3 / (1 - r3)

# The total sum is S = (1/2)*sum1 + (1/2)*sum2 + sum3
term1 = Fraction(1, 2) * sum1
term2 = Fraction(1, 2) * sum2
term3 = sum3

total_sum = term1 + term2 + term3

# Print the final equation with each number.
# Using Fraction type ensures exact arithmetic with rational numbers.
print("The final sum is the result of the following equation:")
print(f"1/2 * {sum1} + 1/2 * ({sum2}) + {sum3} = {term1} + ({term2}) + {term3} = {total_sum}")
