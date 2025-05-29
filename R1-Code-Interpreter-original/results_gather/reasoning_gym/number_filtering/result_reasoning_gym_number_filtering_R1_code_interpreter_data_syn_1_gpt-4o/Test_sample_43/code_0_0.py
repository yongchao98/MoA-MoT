# Original list of numbers as strings
numbers = ['16.960', '-91.6731', '34.691', '81.95', '-79.99', '85.38', '-45.08', '21.333', '-46', '-24.3579']

# Convert strings to floats and filter numbers greater than 9.691
filtered_numbers = [num for num in numbers if float(num) > 9.691]

# Print the filtered list
print(filtered_numbers)