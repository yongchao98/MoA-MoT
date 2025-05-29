# Step 1: Convert hexadecimal to decimal and multiply
num1 = 0xc18a
num2 = 0x8eea
product = num1 * num2

# Step 2: Perform the left shift operation
result = product << 2

# Step 3: Convert the result to hexadecimal
hex_result = hex(result)

# Print the final result
print(hex_result)