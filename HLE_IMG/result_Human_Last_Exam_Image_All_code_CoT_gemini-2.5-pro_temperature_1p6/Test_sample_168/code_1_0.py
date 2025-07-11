import math

# This problem can be solved using combinatorics by finding the number of
# valid sequences of reflections.

# Step 1: Define the number of reflections on each type of mirror.
# There are 2 reflections on G1 and 1 on G3. These are vertical mirrors.
num_G1 = 2
num_G3 = 1
num_vertical_reflections = num_G1 + num_G3

# There is 1 reflection on G2 and 1 on G4. These are horizontal mirrors.
num_G2 = 1
num_G4 = 1
num_horizontal_reflections = num_G2 + num_G4

print(f"The total number of reflections on vertical mirrors (G1, G3) is {num_vertical_reflections}.")
print(f"The total number of reflections on horizontal mirrors (G2, G4) is {num_horizontal_reflections}.")
print("-" * 20)

# Step 2: Determine the number of ways to arrange the reflections in their slots.
# For a path to be possible, reflections must alternate between vertical and horizontal mirrors.
# With 3 vertical and 2 horizontal reflections, the sequence must be V-H-V-H-V.

# The number of ways to arrange the vertical reflections {G1, G1, G3} is 3! / 2!.
perms_vertical = math.factorial(num_vertical_reflections) // math.factorial(num_G1)
num_vertical_reflections_str = str(num_vertical_reflections) + "!"
num_G1_str = str(num_G1) + "!"
print(f"Number of ways to arrange the vertical reflections: {num_vertical_reflections_str} / {num_G1_str} = {perms_vertical}")


# The number of ways to arrange the horizontal reflections {G2, G4} is 2!.
perms_horizontal = math.factorial(num_horizontal_reflections)
num_horizontal_reflections_str = str(num_horizontal_reflections) + "!"
print(f"Number of ways to arrange the horizontal reflections: {num_horizontal_reflections_str} = {perms_horizontal}")
print("-" * 20)

# Step 3: Calculate the total number of unique paths.
# This is the product of the two permutations calculated above.
total_ways = perms_vertical * perms_horizontal

# Print the final calculation and the result.
print("The total number of ways is the product of the possible arrangements for each group.")
print(f"Final Equation: {perms_vertical} * {perms_horizontal} = {total_ways}")
<<<6>>>