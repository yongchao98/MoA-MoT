# Define the numbers as strings to handle commas
numbers = ["312171335.652000", "312171339.704", "312171331.531"]

# Convert the strings to floats for comparison
numbers = [float(num.replace(',', '')) for num in numbers]

# Find the largest number
largest_number = max(numbers)

# Print the largest number
print(largest_number)