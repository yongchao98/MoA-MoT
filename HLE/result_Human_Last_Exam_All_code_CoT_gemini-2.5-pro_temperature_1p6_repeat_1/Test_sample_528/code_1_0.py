# Step 1: Determine the number of choices for placing each element based on the constraints.

# For element 1: must be in S1. Only 1 choice (Region 1).
choices_for_1 = 1

# For element 2: must be in S2. 2 choices (Region 1 or Region 2).
choices_for_2 = 2

# For element 3: must be in S3. 3 choices (Region 1, 2, or 3).
choices_for_3 = 3

# For element 4: no constraints. 5 choices (any of the 5 regions).
choices_for_4 = 5

# For element 5: no constraints. 5 choices (any of the 5 regions).
choices_for_5 = 5

# Step 2: Calculate the total number of sets by multiplying the choices for each element.
total_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

# Step 3: Print the calculation and the final result.
# The final result is the product of the number of choices for each element.
print(f"The number of sets is found by multiplying the number of choices for each element:")
print(f"{choices_for_1} (for 1) * {choices_for_2} (for 2) * {choices_for_3} (for 3) * {choices_for_4} (for 4) * {choices_for_5} (for 5) = {total_sets}")