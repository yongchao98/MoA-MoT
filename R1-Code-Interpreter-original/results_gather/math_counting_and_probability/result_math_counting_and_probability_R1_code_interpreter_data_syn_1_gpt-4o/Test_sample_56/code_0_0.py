# Possible values for r and k
r_values = [-2, -1, 0, 1, 2, 3, 4, 5]
k_values = [2, 3, 4, 5, 6, 7]

# Count favorable outcomes
favorable_outcomes = 0

for k in k_values:
    for r in r_values:
        if r % k == 0:
            favorable_outcomes += 1

# Total outcomes
total_outcomes = len(r_values) * len(k_values)

# Probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)