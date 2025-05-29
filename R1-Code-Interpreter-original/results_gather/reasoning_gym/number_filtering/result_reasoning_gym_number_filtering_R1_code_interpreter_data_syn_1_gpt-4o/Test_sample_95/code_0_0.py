# Original list of strings
original_list = ['82.4730', '13', '-95.29']

# Convert strings to floats and filter out numbers greater than 71.829
filtered_list = [x for x in original_list if float(x) <= 71.829]

# Print the filtered list
print(filtered_list)