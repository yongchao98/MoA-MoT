from fractions import Fraction

# Sum for pairs (d, d), d >= 1
# This is sum_{d=1 to inf} (1/2^(2d)) = sum (1/4)^d
sum1_val = (Fraction(1, 4)) / (1 - Fraction(1, 4))

# Sum for pairs (d, 2d), d >= 1
# This is sum_{d=1 to inf} (1/2^(3d)) = sum (1/8)^d
sum2_val = (Fraction(1, 8)) / (1 - Fraction(1, 8))

# Sum for pairs (2d, d), d >= 1
# This is sum_{d=1 to inf} (1/2^(3d)) = sum (1/8)^d
sum3_val = (Fraction(1, 8)) / (1 - Fraction(1, 8))

total_sum = sum1_val + sum2_val + sum3_val

print(f"The sum over pairs (d,d) is: {sum1_val}")
print(f"The sum over pairs (d,2d) is: {sum2_val}")
print(f"The sum over pairs (2d,d) is: {sum3_val}")
print("\nThe final sum is the sum of these three parts:")
print(f"{sum1_val} + {sum2_val} + {sum3_val} = {total_sum}")