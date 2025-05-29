# Step 1: Convert base-3 to decimal
base_3_number = "220020"
decimal_number = sum(int(digit) * (3 ** idx) for idx, digit in enumerate(reversed(base_3_number)))

# Step 2: Convert decimal to binary
binary_number = bin(decimal_number)[2:]

print(binary_number)