# Initial values
a = 0x23e8
b = 0x386e

# Perform the right shift operation
b_shifted = b >> 0x3

# Perform the multiplication
result = a * b_shifted

# Print the result in hexadecimal format
print(hex(result))