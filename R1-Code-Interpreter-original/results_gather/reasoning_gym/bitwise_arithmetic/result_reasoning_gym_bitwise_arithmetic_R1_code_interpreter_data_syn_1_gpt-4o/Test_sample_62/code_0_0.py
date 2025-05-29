# Define the hexadecimal values
a = 0x6f42
b = 0xda47
c = 0x5ff3
shift_amount = 0x2

# Perform the operations
multiplication_result = a * b
shift_result = c << shift_amount
final_result = multiplication_result + shift_result

# Print the final result in hexadecimal format
print(hex(final_result))