# Step 1: Multiplication
multiplication_result = 0xb7e1 * 0x6f4

# Step 2: Left Shift by 3 bits
shifted_result = multiplication_result << 3

# Step 3: Convert to hexadecimal
hex_result = hex(shifted_result)

# Output the result
print(hex_result)