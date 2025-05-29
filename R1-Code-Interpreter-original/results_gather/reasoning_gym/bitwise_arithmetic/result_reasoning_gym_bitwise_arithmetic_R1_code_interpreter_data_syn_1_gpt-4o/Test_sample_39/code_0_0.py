# Convert hexadecimal to decimal
num1 = 0xfc9c
num2 = 0x28a0

# Perform multiplication
product = num1 * num2

# Right shift the result by 1
result = product >> 1

# Convert the result to hexadecimal
hex_result = hex(result)

# Print the final result
print(hex_result)