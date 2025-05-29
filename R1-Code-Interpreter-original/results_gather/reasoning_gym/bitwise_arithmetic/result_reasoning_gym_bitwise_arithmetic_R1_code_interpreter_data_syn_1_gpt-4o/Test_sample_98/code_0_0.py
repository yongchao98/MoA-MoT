# Define the hexadecimal numbers
a = 0xa6b1
b = 0x57a0

# Perform the bitwise operations
left_shift_result = a << 0x0
right_shift_result = b >> 0x2

# Multiply the results
result = left_shift_result * right_shift_result

# Print the result in hexadecimal format
print(hex(result))