# Initial hexadecimal value
value = 0xad67

# Right shift by 3 bits
right_shifted = value >> 0x3

# Left shift the result by 2 bits
result = right_shifted << 0x2

# Print the final result in hexadecimal format
print(hex(result))