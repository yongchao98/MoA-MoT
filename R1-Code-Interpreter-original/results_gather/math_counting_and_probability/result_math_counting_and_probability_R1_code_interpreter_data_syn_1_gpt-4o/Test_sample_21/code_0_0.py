from math import comb

# Given values
n = 4  # number of days
k = 1  # number of sunny days
p = 0.25  # probability of a sunny day

# Calculate the probability using the binomial formula
probability = comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

# Convert the probability to a fraction
from fractions import Fraction
probability_fraction = Fraction(probability).limit_denominator()

print(probability_fraction)