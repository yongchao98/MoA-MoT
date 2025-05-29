# Step 1: Convert base-3 to decimal
base_3_number = "220020"
decimal_number = 0
base = 3

# Convert base-3 to decimal
for i, digit in enumerate(reversed(base_3_number)):
    decimal_number += int(digit) * (base ** i)

# Step 2: Convert decimal to binary
binary_number = ""
while decimal_number > 0:
    binary_number = str(decimal_number % 2) + binary_number
    decimal_number //= 2

# Print the binary number
print(binary_number)