# Initial hexadecimal value
value = 0x3f4d

# Perform the first right shift by 3 bits
first_shift = value >> 0x3

# Perform the second right shift by 2 bits
second_shift = first_shift >> 0x2

# Print the final result in hexadecimal format
print(hex(second_shift))