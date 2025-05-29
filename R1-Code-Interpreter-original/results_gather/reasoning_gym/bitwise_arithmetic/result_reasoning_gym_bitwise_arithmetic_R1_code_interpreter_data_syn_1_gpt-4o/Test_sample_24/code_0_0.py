# Initial value
value = 0x4e17

# First right shift by 0x3
shifted_once = value >> 0x3

# Second right shift by 0x0 (no change)
shifted_twice = shifted_once >> 0x0

# Print the result in hexadecimal format
print(hex(shifted_twice))