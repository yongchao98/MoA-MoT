# Step 1: Create a list of the given numbers
numbers = [11.2, 49.0, 41.54]

# Step 2: Sort the list in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Step 3: Convert each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Step 4: Print the final list of strings
print(sorted_numbers_as_strings)