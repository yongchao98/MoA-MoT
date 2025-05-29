# Original list of strings
original_list = ['15', '-30', '94.39', '45.6363', '-87.069', '82.027', '-93.09', '-99']

# Convert strings to floats and filter based on the condition
filtered_list = [num for num in original_list if float(num) < 46.60]

# Print the filtered list
print(filtered_list)