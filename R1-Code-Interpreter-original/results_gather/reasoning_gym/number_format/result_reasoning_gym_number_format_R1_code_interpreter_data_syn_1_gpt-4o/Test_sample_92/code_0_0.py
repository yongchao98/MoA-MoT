# List of numbers as strings
numbers = ["548,255,761.25", "548255754.792000", "548,255,770.56", "548255759.574000"]

# Remove commas and convert to float
numbers_float = [float(num.replace(',', '')) for num in numbers]

# Find the smallest number
smallest_number = min(numbers_float)

# Print the smallest number
print(smallest_number)