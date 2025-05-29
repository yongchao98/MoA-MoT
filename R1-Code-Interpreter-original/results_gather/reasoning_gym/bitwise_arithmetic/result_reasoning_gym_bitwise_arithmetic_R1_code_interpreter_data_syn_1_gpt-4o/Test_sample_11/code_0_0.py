# Step 1: Perform the multiplication
multiplication_result = 0xacf1 * 0xb3cc

# Step 2: Perform the left shift
left_shift_result = 0x9a4b << 0x0

# Step 3: Subtract the left shift result from the multiplication result
final_result = multiplication_result - left_shift_result

# Output the final result in hexadecimal format
print(hex(final_result))