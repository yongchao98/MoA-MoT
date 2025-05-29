# Original list of strings
original_list = ['-80.80', '49', '31.103', '54.4', '84.485', '-24.2341']

# Convert strings to floats and filter out numbers larger than -79.2
filtered_list = [x for x in original_list if float(x) <= -79.2]

# Print the filtered list
print(filtered_list)