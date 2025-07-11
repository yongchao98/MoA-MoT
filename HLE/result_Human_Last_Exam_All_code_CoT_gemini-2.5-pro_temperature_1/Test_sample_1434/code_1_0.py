from fractions import Fraction

# The sum can be broken down into three geometric series.
# Sum = 1/2 * Sum_1 + 1/2 * Sum_2 + Sum_3

# Sum_1 = sum_{n=1 to inf} (3/4)^n
r1 = Fraction(3, 4)
sum_1 = r1 / (1 - r1)

# Sum_2 = sum_{n=1 to inf} (-1/4)^n
r2 = Fraction(-1, 4)
sum_2 = r2 / (1 - r2)

# Sum_3 = sum_{n=1 to inf} (1/4)^n
r3 = Fraction(1, 4)
sum_3 = r3 / (1 - r3)

# The total sum is (N_n + 1) / 4^n summed over n, which simplifies to:
# 1/2 * Sum_1 + 1/2 * Sum_2 + Sum_3
total_sum = Fraction(1, 2) * sum_1 + Fraction(1, 2) * sum_2 + sum_3

print("The calculation is based on breaking the sum into three parts:")
print(f"Sum_1 = sum( (3/4)^n ) = {sum_1}")
print(f"Sum_2 = sum( (-1/4)^n ) = {sum_2}")
print(f"Sum_3 = sum( (1/4)^n ) = {sum_3}")
print("\nThe final sum is:")
print(f"(1/2) * Sum_1 + (1/2) * Sum_2 + Sum_3 = (1/2) * {sum_1} + (1/2) * ({sum_2}) + {sum_3} = {total_sum}")