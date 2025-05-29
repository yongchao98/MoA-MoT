# Define the hexadecimal numbers
num1 = 0xff52
num2 = 0x80ac

# Perform the right shift and left shift operations
right_shift_result = num1 >> 0x2
left_shift_result = num2 << 0x2

# Add the results
result = right_shift_result + left_shift_result

# Print the result in hexadecimal format
print(hex(result))