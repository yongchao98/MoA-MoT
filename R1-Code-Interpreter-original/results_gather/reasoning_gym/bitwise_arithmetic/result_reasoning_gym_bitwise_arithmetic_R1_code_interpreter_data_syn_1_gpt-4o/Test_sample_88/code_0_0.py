# Step 1: Right shift 0xacc3 by 3
right_shift_result = 0xacc3 >> 3

# Step 2: Subtract 0x3a20 from 0xdc9b
subtraction_result = 0x3a20 - 0xdc9b

# Step 3: Subtract the result of step 2 from step 1
final_result = right_shift_result - subtraction_result

# Output the final result in hexadecimal format
print(hex(final_result))