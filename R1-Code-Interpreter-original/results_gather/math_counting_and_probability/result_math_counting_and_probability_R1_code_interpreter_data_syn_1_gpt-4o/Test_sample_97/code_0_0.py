# Total number of ways to paint the cube
total_paintings = 2**6

# Favorable outcomes: 2 scenarios for top and bottom, each with 2 options for sides
favorable_outcomes = 2 * 2

# Probability calculation
probability = favorable_outcomes / total_paintings

print(probability)