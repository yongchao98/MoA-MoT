from itertools import product

def is_composite(n):
    if n <= 1:
        return False
    if n <= 3:
        return False
    if n % 2 == 0 or n % 3 == 0:
        return True
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return True
        i += 6
    return False

# Total possible outcomes
total_outcomes = 6 ** 5

# Generate all possible outcomes
dice_faces = [1, 2, 3, 4, 5, 6]
outcomes = product(dice_faces, repeat=5)

# Count composite outcomes
composite_count = 0
for outcome in outcomes:
    product_value = 1
    for number in outcome:
        product_value *= number
    if is_composite(product_value):
        composite_count += 1

# Calculate probability
probability_composite = composite_count / total_outcomes
print(probability_composite)