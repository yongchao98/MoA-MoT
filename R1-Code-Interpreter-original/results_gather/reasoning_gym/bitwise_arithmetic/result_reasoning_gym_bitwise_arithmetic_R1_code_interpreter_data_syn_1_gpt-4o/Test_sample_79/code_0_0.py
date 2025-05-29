# Step 1: Multiplication
multiplication_result = 0x7f71 * 0x32b7

# Step 2: Bitwise Right Shift
right_shift_result = 0xd0c9 >> 0x1

# Step 3: Subtraction
final_result = multiplication_result - right_shift_result

# Convert the result to hexadecimal
hex_result = hex(final_result)

# Print the final result
print(hex_result)