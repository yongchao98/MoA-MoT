# Number of distinct composite transformations for the vertical mirrors (G1, G1, G3)
num_ways_vertical = 2

# Number of distinct composite transformations for the horizontal mirrors (G2, G4)
num_ways_horizontal = 2

# The total number of ways is the product of the two independent possibilities.
total_ways = num_ways_vertical * num_ways_horizontal

# Print the final calculation and the result.
print(f"The number of distinct transformations for the vertical part is: {num_ways_vertical}")
print(f"The number of distinct transformations for the horizontal part is: {num_ways_horizontal}")
print("The total number of ways to draw the light ray path is the product of these two numbers.")
print(f"Total ways = {num_ways_vertical} * {num_ways_horizontal} = {total_ways}")
