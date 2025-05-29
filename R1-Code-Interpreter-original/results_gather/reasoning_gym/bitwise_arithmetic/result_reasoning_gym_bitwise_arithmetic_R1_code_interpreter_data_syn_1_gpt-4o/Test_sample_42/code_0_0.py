# Step 1: Convert hex to decimal and multiply
result = (0x1fe0 * 0xcdcc)

# Step 2: Right shift the result by 3
result >>= 0x3

# Step 3: Convert the result to hexadecimal
hex_result = hex(result)

# Print the final result
print(hex_result)