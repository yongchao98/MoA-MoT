from math import gcd

# Total number of outcomes
total_outcomes = 20 * 17 * 11 * 3

# Number of favorable outcomes (rolling 19 or 20 on 1d20)
favorable_outcomes = 2 * 17 * 11 * 3

# Calculate the probability as a fraction
numerator = favorable_outcomes
denominator = total_outcomes

# Reduce the fraction
common_divisor = gcd(numerator, denominator)
reduced_numerator = numerator // common_divisor
reduced_denominator = denominator // common_divisor

print(f"{reduced_numerator}/{reduced_denominator}")