# Define the hexadecimal numbers
a = 0xb59c
b = 0x659d

# Perform the subtraction
result = a - b

# Perform the left shift by 0 (which does not change the number)
shifted_result = result << 0x0

# Print the result in hexadecimal format
print(hex(shifted_result))