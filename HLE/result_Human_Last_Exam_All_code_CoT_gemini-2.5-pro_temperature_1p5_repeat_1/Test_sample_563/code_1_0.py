# Define the number of automorphism groups for each genus as specified in the problem.
num_groups_genus_2 = 12
num_groups_genus_3 = 36
num_groups_genus_4 = 23

# Create a list containing these numbers.
result_list = [num_groups_genus_2, num_groups_genus_3, num_groups_genus_4]

# Print the list in the specified format.
# The f-string explicitly uses each number to construct the output string.
print(f"[{result_list[0]},{result_list[1]},{result_list[2]}]")