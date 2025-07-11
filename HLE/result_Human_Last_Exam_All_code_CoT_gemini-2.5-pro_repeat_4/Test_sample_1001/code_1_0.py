from fractions import Fraction

def geometric_series_sum(first_term, ratio):
    """
    Calculates the sum of an infinite geometric series.
    The formula is S = a / (1 - r), where a is the first term.
    """
    if abs(ratio) >= 1:
        raise ValueError("The absolute value of the ratio must be less than 1 for the series to converge.")
    return first_term / (1 - ratio)

# The total sum is the sum of three distinct series, corresponding to the three
# families of pairs (i, j) in the set S.

# 1. For pairs of the form (d, d), the sum is Sum_{d=1 to inf} 1/2^(2d).
# This is a geometric series with first term a = 1/4 and ratio r = 1/4.
sum1 = geometric_series_sum(Fraction(1, 4), Fraction(1, 4))

# 2. For pairs of the form (d, 2d), the sum is Sum_{d=1 to inf} 1/2^(3d).
# This is a geometric series with first term a = 1/8 and ratio r = 1/8.
sum2 = geometric_series_sum(Fraction(1, 8), Fraction(1, 8))

# 3. For pairs of the form (2d, d), the sum is Sum_{d=1 to inf} 1/2^(3d).
# This is the same series as for the (d, 2d) case.
sum3 = geometric_series_sum(Fraction(1, 8), Fraction(1, 8))

# The total sum is the sum of these three parts.
total_sum = sum1 + sum2 + sum3

# Print the final equation with each component clearly stated.
print(f"The required sum is composed of three parts:")
print(f"Sum for pairs (d, d) = {sum1}")
print(f"Sum for pairs (d, 2d) = {sum2}")
print(f"Sum for pairs (2d, d) = {sum3}")
print("\nThe final equation is:")
print(f"{sum1} + {sum2} + {sum3} = {total_sum}")
