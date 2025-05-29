# Initial value
value = 0x4550

# First shift: left shift by 3
result_after_first_shift = value << 0x3

# Second shift: left shift by 1
final_result = result_after_first_shift << 0x1

# Print the final result in hexadecimal format
print(hex(final_result))