# Step 1: Perform the subtraction
result_subtraction = 0x6646 - 0xbb94

# Step 2: Perform the left shift
result_shift = result_subtraction << 1

# Step 3: Convert the result to hexadecimal
# If the result is negative, it should be represented as a negative hexadecimal value
if result_shift < 0:
    result_hex = f"-0x{-result_shift:x}"
else:
    result_hex = f"0x{result_shift:x}"

print(result_hex)