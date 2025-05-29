from itertools import combinations

# Define the set
numbers = [2, 4, 10, 12, 15, 20, 50]

# Function to check if a product is a multiple of 100
def is_multiple_of_100(a, b):
    product = a * b
    return product % 100 == 0

# Calculate total number of combinations
total_combinations = len(list(combinations(numbers, 2)))

# Calculate number of favorable combinations
favorable_combinations = sum(1 for a, b in combinations(numbers, 2) if is_multiple_of_100(a, b))

# Calculate probability
probability = favorable_combinations / total_combinations

print(probability)