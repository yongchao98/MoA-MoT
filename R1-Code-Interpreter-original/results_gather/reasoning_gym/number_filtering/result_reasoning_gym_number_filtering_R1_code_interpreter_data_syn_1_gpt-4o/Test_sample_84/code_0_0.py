# Original list of numbers as strings
numbers = ['-98', '47.525', '93.5', '-45.1788', '-15.2636', '-85', '-36.9657']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) <= -36]

# Print the new list
print(filtered_numbers)