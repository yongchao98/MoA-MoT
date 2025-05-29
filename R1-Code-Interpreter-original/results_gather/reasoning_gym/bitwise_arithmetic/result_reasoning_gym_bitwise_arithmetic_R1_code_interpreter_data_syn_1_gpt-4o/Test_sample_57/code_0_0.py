# Step 1: Multiply the two hexadecimal numbers
result = 0x24fd * 0x5922

# Step 2: Right shift the result by 2 bits (equivalent to dividing by 4)
result >>= 2

# Step 3: Convert the result to hexadecimal
# If the result is negative, Python will automatically handle the sign
hex_result = hex(result)

# Print the final result
print(hex_result)