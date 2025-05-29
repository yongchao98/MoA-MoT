# Original list of strings
original_list = ['11.24', '49', '41.54']

# Convert strings to floats and filter
filtered_list = [num for num in original_list if float(num) > 42.255]

# Print the filtered list
print(filtered_list)