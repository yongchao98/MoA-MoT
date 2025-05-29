from itertools import product

# Possible numbers on a die
dice_numbers = [1, 2, 3, 4, 5, 6]

# Count the number of favorable outcomes
favorable_count = 0

# Iterate over all possible outcomes of three dice
for a, b, c in product(dice_numbers, repeat=3):
    if a * b * c == 180:
        favorable_count += 1

# Total possible outcomes
total_outcomes = 6 ** 3

# Calculate the probability
probability = favorable_count / total_outcomes

print(probability)