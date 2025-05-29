# Original list of strings
original_list = ['-62.2704', '34.189', '-81.8', '-55.37', '44.58', '48']

# Convert strings to floats and filter numbers smaller than 9.87
filtered_list = [num for num in original_list if float(num) < 9.87]

# Print the filtered list
print(filtered_list)