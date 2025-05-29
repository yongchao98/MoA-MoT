from itertools import product

# Possible outcomes for a single die
outcomes = [1, 2, 3, 4, 5, 6]

# Count the number of favorable outcomes
favorable_count = 0
total_count = 0

# Iterate over all possible outcomes for four dice
for dice in product(outcomes, repeat=4):
    product_even = any(d % 2 == 0 for d in dice)  # Check if product is even
    sum_odd = sum(dice) % 2 == 1  # Check if sum is odd
    
    if product_even:
        total_count += 1
        if sum_odd:
            favorable_count += 1

# Calculate the probability
probability = favorable_count / total_count
print(probability)