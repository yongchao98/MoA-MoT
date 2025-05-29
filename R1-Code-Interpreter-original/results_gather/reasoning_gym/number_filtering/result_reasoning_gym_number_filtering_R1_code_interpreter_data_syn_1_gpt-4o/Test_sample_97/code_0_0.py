# Original list of strings
original_list = ['-90.765', '26.17', '73.1', '-77.6256', '-69.37', '13.2', '-3.7093']

# Convert strings to floats and filter out numbers greater than 39.4
filtered_list = [num for num in original_list if float(num) <= 39.4]

# Print the filtered list
print(filtered_list)