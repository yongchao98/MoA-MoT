# Original list of strings
original_list = ['-5.2656', '-39.7718', '97.721']

# Convert strings to floats and filter numbers smaller than 57
filtered_list = [num_str for num_str in original_list if float(num_str) < 57]

# Print the filtered list
print(filtered_list)