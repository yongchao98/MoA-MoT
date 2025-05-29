# Define the numbers as strings to handle commas
numbers = ["737915332.359000", "737915339.862", "737915338.624"]

# Remove commas and convert to float
numbers = [float(num.replace(',', '')) for num in numbers]

# Find the smallest number
smallest_number = min(numbers)

# Print the smallest number
print(smallest_number)