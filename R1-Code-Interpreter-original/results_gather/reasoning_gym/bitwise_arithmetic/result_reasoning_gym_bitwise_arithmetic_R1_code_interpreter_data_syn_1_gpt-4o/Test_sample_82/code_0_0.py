# Step 1: Subtraction
subtraction_result = 0xa154 - 0xa0cf

# Step 2: Right Shift
right_shift_result = 0xfc16 >> 0x2

# Step 3: Addition
final_result = subtraction_result + right_shift_result

# Convert the final result to hexadecimal
hex_result = hex(final_result)

# Print the final result
print(hex_result)