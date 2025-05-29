# Step 1: Subtract 0x533b from 0x849d
result_subtraction = 0x849d - 0x533b

# Step 2: Right shift the result by 3 bits
result_shift = result_subtraction >> 3

# Print the final result in hexadecimal format
print(hex(result_shift))