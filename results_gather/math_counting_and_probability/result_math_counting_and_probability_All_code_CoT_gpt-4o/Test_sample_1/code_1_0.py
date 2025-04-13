from math import factorial

# Total number of arrangements without considering symmetry
total_arrangements = factorial(10)

# Number of symmetries (5 rotations + 5 reflections)
num_symmetries = 10

# Using Burnside's Lemma
distinct_arrangements = total_arrangements / num_symmetries

# Output the result
print(f"<<<{distinct_arrangements}>>>")