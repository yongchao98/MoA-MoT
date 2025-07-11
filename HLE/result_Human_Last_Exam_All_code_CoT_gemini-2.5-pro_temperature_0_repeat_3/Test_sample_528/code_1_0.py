# Define the number of choices for each element based on the problem's constraints.
choices_for_1 = 1
choices_for_2 = 2
choices_for_3 = 3
choices_for_4 = 5
choices_for_5 = 5

# Calculate the total number of sets by multiplying the choices for each element.
total_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

# Print the explanation and the final calculation.
print("The number of possible sets (S1, S2, S3, S4) is found by multiplying the number of placement choices for each element in the universal set {1, 2, 3, 4, 5}.")
print(f"Choices for element 1 (must be in S1): {choices_for_1}")
print(f"Choices for element 2 (must be in S2): {choices_for_2}")
print(f"Choices for element 3 (must be in S3): {choices_for_3}")
print(f"Choices for element 4 (no constraint): {choices_for_4}")
print(f"Choices for element 5 (no constraint): {choices_for_5}")
print("\nThe total number of sets is the product of these choices:")
print(f"{choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_sets}")
