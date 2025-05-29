# Define the initial hexadecimal values
value1 = 0x4e0a
value2 = 0x4f57

# Perform the right shift operations
shifted_value1 = value1 >> 1
shifted_value2 = value2 >> 0

# Calculate the result of the subtraction
result = shifted_value1 - shifted_value2

# Print the result in hexadecimal format
# If the result is negative, it will be displayed as a negative value
print(hex(result))