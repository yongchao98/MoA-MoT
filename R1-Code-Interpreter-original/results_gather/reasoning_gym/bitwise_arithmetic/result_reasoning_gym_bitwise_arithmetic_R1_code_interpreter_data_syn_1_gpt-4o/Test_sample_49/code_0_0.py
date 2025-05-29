# Initial value
value = 0x635c

# Right shift by 0 (no change)
shifted_value = value >> 0x0

# Left shift by 2 (equivalent to multiplying by 4)
result = shifted_value << 0x2

# Print the result in hexadecimal format
print(hex(result))