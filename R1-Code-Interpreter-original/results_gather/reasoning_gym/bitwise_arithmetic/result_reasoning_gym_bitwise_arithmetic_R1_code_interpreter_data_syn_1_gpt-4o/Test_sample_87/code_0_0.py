# Define the hexadecimal numbers
num1 = 0xab53
num2 = 0xa273

# Perform the bitwise operations
result1 = num1 >> 0x3
result2 = num2 << 0x1

# Calculate the final result
final_result = result1 - result2

# Print the result in hexadecimal format
print(hex(final_result))