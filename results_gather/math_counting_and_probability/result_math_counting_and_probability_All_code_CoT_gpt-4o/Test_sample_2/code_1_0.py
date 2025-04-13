from math import comb

# Total number of shoes
total_shoes = 22

# Total ways to pick any two shoes
total_ways = comb(total_shoes, 2)

# Ways to pick two shoes of the same color and different sides
black_pairs = 6
brown_pairs = 3
gray_pairs = 2

same_color_different_sides = black_pairs + brown_pairs + gray_pairs

# Probability calculation
probability = same_color_different_sides / total_ways

# Output the result
print(f"<<<{probability}>>>")