# Initial hexadecimal number
number = 0xb766

# First right shift by 1
first_shift = number >> 0x1

# Second right shift by 1
second_shift = first_shift >> 0x1

# Output the result in hexadecimal format
print(hex(second_shift))