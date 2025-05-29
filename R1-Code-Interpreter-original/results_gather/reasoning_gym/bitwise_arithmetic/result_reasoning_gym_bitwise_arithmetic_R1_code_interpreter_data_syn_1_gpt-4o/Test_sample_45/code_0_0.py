# Step 1: Subtraction
subtraction_result = 0xd393 - 0xf18d

# Step 2: Left Shift
left_shift_result = 0xf4fa << 1

# Step 3: Multiplication
multiplication_result = subtraction_result * left_shift_result

# Convert the final result to hexadecimal
final_result_hex = hex(multiplication_result)

# Print the final result
print(final_result_hex)