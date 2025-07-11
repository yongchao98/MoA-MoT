# Plan:
# 1. Define the number of choices for each segment of the path from A to B.
#    The path is conceptually broken into three parts: A->C, C->D, and D->B.
# 2. The number of ways to travel from A to C is 1 (along the line segment).
# 3. The number of ways to travel from C to D is 3 (along the line segment, or along one of the two circle arcs).
# 4. The number of ways to travel from D to B is 1 (along the line segment).
# 5. The total number of distinct paths is the product of these choices.
# 6. The program will print the equation and the final result.

# Number of paths for each segment
choices_A_to_C = 1
choices_C_to_D = 3
choices_D_to_B = 1

# Calculate the total number of paths by multiplying the choices
total_paths = choices_A_to_C * choices_C_to_D * choices_D_to_B

# Print the final equation and the answer
print(f"The number of distinct paths is found by multiplying the choices at each stage:")
print(f"{choices_A_to_C} (A->C) * {choices_C_to_D} (C->D) * {choices_D_to_B} (D->B) = {total_paths}")

print(f"\nTherefore, there are {total_paths} distinct paths.")