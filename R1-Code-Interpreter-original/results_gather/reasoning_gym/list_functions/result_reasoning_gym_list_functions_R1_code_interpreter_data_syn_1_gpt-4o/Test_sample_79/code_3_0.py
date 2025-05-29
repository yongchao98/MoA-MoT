# Define the input list
input_list = [41, 87, 98]

# Sort the list in descending order
sorted_desc = sorted(input_list, reverse=True)

# Sort the list in ascending order
sorted_asc = sorted(input_list)

# Determine the output list by combining the two sorted lists
# Since the input list has only three elements, we can directly use the sorted lists
output_list = sorted_desc[:2] + sorted_asc[2:]

# Print the output list
print(output_list)