from fractions import Fraction

# Step 1: Define the sum based on the characterization of the set S.
# Based on mathematical analysis, for each integer n >= 1, the set of
# possible values for m contains exactly two elements: 0 and N_n, where
# N_n is the size of the neighborhood S_p.
#
# This means the sum can be written as:
# Sum = Sum_{n=1 to infinity} (for m=0) 1/4^n + (for m=N_n) 1/4^n
# Sum = Sum_{n=1 to infinity} 2 * (1/4)^n
# Sum = 2 * Sum_{n=1 to infinity} (1/4)^n

# Step 2: Calculate the value of the geometric series.
# The series is Sum_{k=1 to infinity} r^k, where r = 1/4.
# The sum of such a series is given by the formula: (first term) / (1 - ratio).
# The first term is 1/4, and the ratio is 1/4.

first_term = Fraction(1, 4)
ratio = Fraction(1, 4)
series_sum = first_term / (1 - ratio)

# The total sum is twice the value of this geometric series.
total_sum = 2 * series_sum

print("The problem asks for the sum of 1/4^n over a set S of pairs (n, m).")
print("For each n >= 1, there are two possible values of m for which a valid set P exists.")
print("Therefore, the required sum is equivalent to:")
print("Sum = 2 * (1/4 + 1/16 + 1/64 + ...)")
print("\nThis is twice a geometric series.")
print(f"The first term of the series is a = {first_term.numerator}/{first_term.denominator}")
print(f"The common ratio is r = {ratio.numerator}/{ratio.denominator}")
print("The sum of the series is a / (1 - r).")

# Output each number in the final equation
print("\nCalculation:")
print(f"Sum of series = ({first_term.numerator}/{first_term.denominator}) / (1 - {ratio.numerator}/{ratio.denominator})")
interim_denominator = 1 - ratio
print(f"             = ({first_term.numerator}/{first_term.denominator}) / ({interim_denominator.numerator}/{interim_denominator.denominator})")
print(f"             = {series_sum.numerator}/{series_sum.denominator}")

print("\nThe final sum is twice the sum of the series:")
print(f"Final Sum = 2 * {series_sum.numerator}/{series_sum.denominator} = {total_sum.numerator}/{total_sum.denominator}")
print(f"\nThe result is {total_sum}")
