# Step 1: Convert base-7 to decimal
base_7_number = "1051"
decimal_number = 0
base = 7

# Convert base-7 to decimal
for i, digit in enumerate(reversed(base_7_number)):
    decimal_number += int(digit) * (base ** i)

# Step 2: Convert decimal to hexadecimal
hexadecimal_number = ""
hex_base = 16
digits = "0123456789abcdef"

while decimal_number > 0:
    remainder = decimal_number % hex_base
    hexadecimal_number = digits[remainder] + hexadecimal_number
    decimal_number //= hex_base

print(hexadecimal_number)