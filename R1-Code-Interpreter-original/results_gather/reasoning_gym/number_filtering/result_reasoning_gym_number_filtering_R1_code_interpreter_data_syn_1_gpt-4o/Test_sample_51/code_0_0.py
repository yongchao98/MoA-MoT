# Original list of numbers as strings
numbers = ['68.759', '-82.261', '-93.1', '-62.16', '-83', '-17.87', '-39.7', '-4.9', '64.3', '-91']

# Convert strings to floats and filter out numbers smaller than -89.886
filtered_numbers = [num for num in numbers if float(num) >= -89.886]

# Print the filtered list
print(filtered_numbers)