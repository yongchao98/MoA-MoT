# Step 1: Subtraction
result_subtraction = 0xd33b - 0xff5d

# Step 2: Right Shift
result_shift = result_subtraction >> 0x3

# Output the final result in hexadecimal
print(hex(result_shift))