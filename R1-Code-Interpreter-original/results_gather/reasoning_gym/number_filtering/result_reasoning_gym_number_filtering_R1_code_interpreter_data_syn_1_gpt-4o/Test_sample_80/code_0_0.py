# Original list of numbers as strings
numbers = ['66.5934', '23.967', '80.5', '-27.5017', '-9.041', '19.1', '45.744', '-47.5', '-54.5255', '34.7']

# Convert strings to floats and filter out numbers smaller than 46.2790
filtered_numbers = [num for num in numbers if float(num) >= 46.2790]

# Print the filtered list
print(filtered_numbers)