# Initial value
value = 0x58ed

# Right shift by 3
right_shifted = value >> 0x3

# Left shift by 0 (no change)
result = right_shifted << 0x0

# Print the result in hexadecimal format
print(hex(result))