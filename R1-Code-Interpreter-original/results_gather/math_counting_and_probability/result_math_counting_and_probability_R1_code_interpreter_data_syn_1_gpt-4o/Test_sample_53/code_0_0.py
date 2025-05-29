# Total number of balls
total_balls = 50

# Count multiples of 5
multiples_of_5 = len([i for i in range(1, total_balls + 1) if i % 5 == 0])

# Count multiples of 7
multiples_of_7 = len([i for i in range(1, total_balls + 1) if i % 7 == 0])

# Count multiples of both 5 and 7 (i.e., multiples of 35)
multiples_of_35 = len([i for i in range(1, total_balls + 1) if i % 35 == 0])

# Apply the principle of inclusion-exclusion
probability = (multiples_of_5 + multiples_of_7 - multiples_of_35) / total_balls

# Print the probability as a fraction
from fractions import Fraction
print(Fraction(probability).limit_denominator())