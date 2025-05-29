# Define the hexadecimal values
right_shift_result = 0x4e60 >> 0x1
left_shift_result = 0x2b67 << 0x0

# Perform the subtraction
result = right_shift_result - left_shift_result

# Print the result in hexadecimal format
print(hex(result))