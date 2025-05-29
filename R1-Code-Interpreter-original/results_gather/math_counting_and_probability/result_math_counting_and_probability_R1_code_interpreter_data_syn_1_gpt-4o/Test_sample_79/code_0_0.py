from itertools import combinations

# Stick lengths
sticks = [2, 3, 5, 7, 11, 13, 17]

# Function to check if three lengths can form a triangle
def can_form_triangle(a, b, c):
    return a + b > c and a + c > b and b + c > a

# Generate all combinations of 3 sticks
all_combinations = combinations(sticks, 3)

# Count the number of valid triangles
valid_triangles = sum(1 for comb in all_combinations if can_form_triangle(*sorted(comb)))

# Total number of combinations
total_combinations = len(list(combinations(sticks, 3)))

# Calculate probability
probability = valid_triangles / total_combinations

print(probability)