import math

# Total number of shoes
total_shoes = 22

# Total ways to pick 2 shoes from 22
total_ways = math.comb(total_shoes, 2)

# Favorable outcomes
favorable_black = 6  # 6 pairs of black shoes
favorable_brown = 3  # 3 pairs of brown shoes
favorable_gray = 2   # 2 pairs of gray shoes

# Total favorable outcomes
favorable_outcomes = favorable_black + favorable_brown + favorable_gray

# Probability
probability = favorable_outcomes / total_ways

# Output the probability as a common fraction
from fractions import Fraction
probability_fraction = Fraction(probability).limit_denominator()

print(probability_fraction)