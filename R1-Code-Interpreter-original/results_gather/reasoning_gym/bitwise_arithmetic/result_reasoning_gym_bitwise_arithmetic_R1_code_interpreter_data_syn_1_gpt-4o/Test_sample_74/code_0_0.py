# Initial value
value = 0xd32d

# Right shift by 0 (no change)
value = value >> 0x0

# Left shift by 1 (equivalent to multiplying by 2)
result = value << 0x1

# Print the result in hexadecimal format
print(hex(result))