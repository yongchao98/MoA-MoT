# Step 1: Subtract 0x6f85 from 0x4a31
step1_result = 0x4a31 - 0x6f85

# Step 2: Right shift 0x10d2 by 2 bits
step2_result = 0x10d2 >> 2

# Step 3: Subtract the result of Step 2 from the result of Step 1
final_result = step1_result - step2_result

# Print the final result in hexadecimal format
print(hex(final_result))