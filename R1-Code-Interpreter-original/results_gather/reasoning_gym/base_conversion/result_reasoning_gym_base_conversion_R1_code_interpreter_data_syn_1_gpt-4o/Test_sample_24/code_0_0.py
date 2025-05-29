# Define the base-15 number
base_15_number = "45e"

# Define a dictionary to map base-15 characters to their integer values
base_15_to_int = {
    '0': 0, '1': 1, '2': 2, '3': 3, '4': 4,
    '5': 5, '6': 6, '7': 7, '8': 8, '9': 9,
    'a': 10, 'b': 11, 'c': 12, 'd': 13, 'e': 14
}

# Initialize the base-10 value
base_10_value = 0

# Calculate the base-10 value
for i, digit in enumerate(reversed(base_15_number)):
    base_10_value += base_15_to_int[digit] * (15 ** i)

# Print the result
print(base_10_value)