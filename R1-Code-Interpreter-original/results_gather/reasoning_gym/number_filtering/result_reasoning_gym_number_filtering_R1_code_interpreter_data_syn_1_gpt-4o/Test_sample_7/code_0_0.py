# Original list of strings
original_list = ['-48', '-60.988', '46', '29.16', '35.3', '92']

# Convert strings to floats and filter numbers smaller than -16
filtered_list = [num_str for num_str in original_list if float(num_str) < -16]

# Print the filtered list
print(filtered_list)