# Define the hexadecimal numbers
num1 = 0xe315
num2 = 0x1e15

# Perform the multiplication
result = num1 * num2

# Perform the left shift (which is effectively no shift here)
shifted_result = result << 0x0

# Convert the result to hexadecimal
hex_result = hex(shifted_result)

# Print the final result
print(hex_result)