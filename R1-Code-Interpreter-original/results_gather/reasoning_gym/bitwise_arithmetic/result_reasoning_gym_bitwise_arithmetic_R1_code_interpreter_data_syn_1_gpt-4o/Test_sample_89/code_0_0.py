# Step 1: Convert hexadecimal to decimal and perform subtraction
result = 0x9e4e - 0xcd8f

# Step 2: Right shift the result by 3 bits
shifted_result = result >> 3

# Step 3: Convert the result to hexadecimal
hex_result = hex(shifted_result)

# Print the final result
print(hex_result)