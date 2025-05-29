# Step 1: Subtract 0x691a from 0x8ffd
result1 = 0x8ffd - 0x691a

# Step 2: Right shift 0xf3a3 by 2 bits
shifted_value = 0xf3a3 >> 2

# Step 3: Subtract the shifted value from the result of step 1
final_result = result1 - shifted_value

# Print the final result in hexadecimal format
print(hex(final_result))