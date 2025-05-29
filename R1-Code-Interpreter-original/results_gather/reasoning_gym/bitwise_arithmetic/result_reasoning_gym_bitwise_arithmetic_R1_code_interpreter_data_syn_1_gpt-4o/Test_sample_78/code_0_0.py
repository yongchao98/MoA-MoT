# Initial value
value = 0x406d

# Right shift by 3
right_shifted = value >> 0x3

# Left shift the result by 2
result = right_shifted << 0x2

# Print the result in hexadecimal format
print(hex(result))