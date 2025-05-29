# Original list of strings
original_list = ['-42.8', '91.88', '34']

# Convert strings to floats and filter out numbers greater than 18.236
filtered_list = [x for x in original_list if float(x) <= 18.236]

# Print the filtered list
print(filtered_list)