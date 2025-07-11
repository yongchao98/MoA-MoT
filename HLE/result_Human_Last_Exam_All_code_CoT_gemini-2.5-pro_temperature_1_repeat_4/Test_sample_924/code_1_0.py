# Number of green cubes in the top layer
# The top face has 3 rows, and each row must have 2 green cubes.
# So, the total number of green cubes in the top layer is 3 * 2.
G_top = 3 * 2

# Number of green cubes in the bottom layer
# Similarly, the bottom layer has 3 rows on its face, each with 2 green cubes.
G_bottom = 3 * 2

# The total number of green cubes is the sum of green cubes in the top, middle, and bottom layers.
# G_total = G_top + G_middle + G_bottom
# G_total = 6 + G_middle + 6 = 12 + G_middle
# We need to find the minimum and maximum possible number of green cubes in the middle layer (G_middle).

# Analysis of the middle layer constraints reveals the possible range for G_middle.
# Let M(x,z) be the color of the cube at (x,1,z), where 1=Green, 0=Red.
# The rules on the side faces provide 4 linear equations for the 9 variables M(x,z).
# Solving this system shows that the sum of all M(x,z), which is G_middle, can be calculated.
# G_middle = 8 + M(1,1) - (M(0,0) + M(0,2) + M(2,0) + M(2,2))
# where M(1,1) is the core cube and the others are the four corner cubes of the middle layer.

# To find the minimum G_middle, we minimize the core cube's value (0) and maximize the corners' sum (4).
G_middle_min = 8 + 0 - 4

# To find the maximum G_middle, we maximize the core cube's value (1) and minimize the corners' sum (2).
G_middle_max = 8 + 1 - 2

# Calculate the smallest and largest total number of green cubes
smallest_total_green = G_top + G_bottom + G_middle_min
largest_total_green = G_top + G_bottom + G_middle_max

print(f"The number of green cubes in the top layer is {G_top}.")
print(f"The number of green cubes in the bottom layer is {G_bottom}.")
print(f"The smallest possible number of green cubes in the middle layer is {G_middle_min}.")
print(f"The largest possible number of green cubes in the middle layer is {G_middle_max}.")
print("-" * 20)
print(f"Smallest total number of green cubes = {G_top} + {G_bottom} + {G_middle_min} = {smallest_total_green}")
print(f"Largest total number of green cubes = {G_top} + {G_bottom} + {G_middle_max} = {largest_total_green}")