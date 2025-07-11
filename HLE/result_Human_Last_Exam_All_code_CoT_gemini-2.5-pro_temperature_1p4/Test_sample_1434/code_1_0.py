from fractions import Fraction

# For n=1, c_1 = 2.
# For n>=2, c_n = (3^n + (-1)^n)/2 + 1.
# The sum is Sum_{n=1 to inf} c_n / 4^n.
# We calculated this to be 26/15. Let's express the calculation step-by-step.

# Contribution for n=1
c1 = 2
term1 = Fraction(c1, 4**1)
print(f"For n=1, the number of possible values for m is c_1 = 2.")
print(f"The term for n=1 in the sum is 2/4^1 = {term1}.")
print("")

# Sum for n>=2
# Sum_{n=2 to inf} ( (3^n+(-1)^n)/2 + 1 ) / 4^n
# = 1/2 * Sum (3/4)^n + 1/2 * Sum (-1/4)^n + Sum (1/4)^n for n from 2 to inf

# Geometric series sum formula for |a|<1: Sum_{n=k to inf} a^n = a^k / (1-a)
# For k=2
def geometric_sum_from_2(a):
    return a**2 / (Fraction(1) - a)

# First part of the sum for n>=2
a1 = Fraction(3, 4)
sum1 = Fraction(1, 2) * geometric_sum_from_2(a1)
print("For n>=2, the sum splits into three geometric series.")
print("Part 1: 1/2 * Sum_{n=2 to inf} (3/4)^n")
print(f"The sum of the geometric series (3/4)^n from n=2 is (3/4)^2 / (1 - 3/4) = {geometric_sum_from_2(a1)}")
print(f"So Part 1 evaluates to 1/2 * {geometric_sum_from_2(a1)} = {sum1}")
print("")

# Second part
a2 = Fraction(-1, 4)
sum2 = Fraction(1, 2) * geometric_sum_from_2(a2)
print("Part 2: 1/2 * Sum_{n=2 to inf} (-1/4)^n")
print(f"The sum of the geometric series (-1/4)^n from n=2 is (-1/4)^2 / (1 - (-1/4)) = {geometric_sum_from_2(a2)}")
print(f"So Part 2 evaluates to 1/2 * {geometric_sum_from_2(a2)} = {sum2}")
print("")

# Third part
a3 = Fraction(1, 4)
sum3 = geometric_sum_from_2(a3)
print("Part 3: Sum_{n=2 to inf} (1/4)^n")
print(f"The sum of the geometric series (1/4)^n from n=2 is (1/4)^2 / (1 - 1/4) = {sum3}")
print("")

# Total sum
total_sum = term1 + sum1 + sum2 + sum3
print("The total sum is the sum of the n=1 term and the three parts for n>=2.")
print(f"Total Sum = {term1} + {sum1} + {sum2} + {sum3}")
print(f"Total Sum = {term1 + sum1} + {sum2} + {sum3}")
print(f"Total Sum = {term1 + sum1 + sum2} + {sum3}")
print(f"Total Sum = {total_sum}")
print(f"The value of the sum is {float(total_sum)}.")
print(f"Final Answer in fraction form: {total_sum.numerator}/{total_sum.denominator}")