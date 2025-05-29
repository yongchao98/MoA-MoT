from math import comb

# Total number of outcomes for 6 coin tosses
total_outcomes = 2**6

# Number of favorable outcomes (exactly 2 heads)
favorable_outcomes = comb(6, 2)

# Probability as a fraction
probability = favorable_outcomes / total_outcomes

# Expressing the probability as a common fraction
from fractions import Fraction
probability_fraction = Fraction(probability).limit_denominator()

print(probability_fraction)