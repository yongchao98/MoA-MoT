# Original list of strings
original_list = ['-89.0955', '-3.23', '-18.6', '-68.170', '66.11', '-12', '-93.161']

# Convert strings to floats and filter numbers greater than -23.4
filtered_list = [num for num in original_list if float(num) > -23.4]

# Print the filtered list
print(filtered_list)