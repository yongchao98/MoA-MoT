# Initial value
value = 0xdd1a

# Right shift by 1
shifted_value = value >> 0x1

# Left shift by 0 (no change)
final_value = shifted_value << 0x0

# Convert to hexadecimal
hex_result = hex(final_value)

# Print the result
print(hex_result)