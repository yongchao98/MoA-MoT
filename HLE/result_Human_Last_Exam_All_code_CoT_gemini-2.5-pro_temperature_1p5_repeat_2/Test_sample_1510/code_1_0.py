# A list of the wallpaper groups that have a point group of order 2.
# These correspond to point groups C2 (for p2) and D1 (for pm, pg, cm).
groups_with_point_group_order_2 = ["p2", "pm", "pg", "cm"]

# Calculate the total count
total_count = len(groups_with_point_group_order_2)

# Print the names of the identified groups
print("The wallpaper groups with a point group of order 2 are:")
for group_name in groups_with_point_group_order_2:
    print(f"- {group_name}")

# Create and print the equation representing the count
# This fulfills the requirement to "output each number in the final equation"
equation_parts = ["1"] * total_count
equation_string = " + ".join(equation_parts)

print(f"\nThe total number is the sum of these groups:")
print(f"{equation_string} = {total_count}")
