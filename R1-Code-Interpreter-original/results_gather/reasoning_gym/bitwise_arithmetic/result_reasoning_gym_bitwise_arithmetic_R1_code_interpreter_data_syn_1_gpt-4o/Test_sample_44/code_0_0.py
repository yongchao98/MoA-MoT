# Step 1: Convert hexadecimal to decimal and multiply
result = 0xe82 * 0x6c09

# Step 2: Perform the left shift operation
result <<= 3

# Step 3: Convert the result to hexadecimal
hex_result = hex(result)

# Print the final result
print(hex_result)