# Define the hexadecimal values
value1 = 0xbb70
shift_amount = 0x3
value2 = 0x4da8
value3 = 0x37fe

# Perform the left shift operation
shifted_value = value1 << shift_amount

# Calculate the sum of value2 and value3
sum_value2_3 = value2 + value3

# Calculate the final result
final_result = shifted_value + sum_value2_3

# Print the final result in hexadecimal format
print(hex(final_result))