# The bins are labeled with integers. The marble starts at bin 0.
# The absorbing states are bin 2024 (melt) and bin 2025 (escape).
# We want to find the probability of escaping given it eventually melts or escapes.

# Let n be the current position of the marble, which is very far from the
# absorbing states (n=0 is far from 2024 and 2025).
# The final fate of the marble is determined by a large jump.
# The probability of jumping from n to a target k is proportional to (1/3)^|k-n|.

# Let P_escape be the probability of jumping to 2025 and P_melt be the probability of jumping to 2024.
# The ratio of these probabilities is:
# Ratio = P_melt / P_escape = (1/3)^(2024-n) / (1/3)^(2025-n)
# The (1/3)^(-n) terms cancel out.
# Ratio = (1/3)^2024 / (1/3)^2025 = 3

# Let p be the probability of escape.
# The probability of melting is 3*p.
# Since the marble must either escape or melt, the probabilities sum to 1.
# p + 3*p = 1
# 4*p = 1
# p = 1/4

# We can write this as 1 / (1 + 3).
numerator = 1
denominator_part_1 = 1
denominator_part_2 = 3

# Calculate the probability
probability = numerator / (denominator_part_1 + denominator_part_2)

# We need to output each number in the final equation
print(f"The probability of the marble escaping is given by the expression: {numerator} / ({denominator_part_1} + {denominator_part_2})")
print(f"The calculated probability is: {probability}")
