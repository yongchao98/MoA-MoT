# Step 1: Convert base-5 number to decimal
base_5_number = "11342"
decimal_number = 0
base = 5

# Convert base-5 to decimal
for i, digit in enumerate(reversed(base_5_number)):
    decimal_number += int(digit) * (base ** i)

# Step 2: Convert decimal to hexadecimal
hexadecimal_number = ""
hex_base = 16
hex_digits = "0123456789abcdef"

while decimal_number > 0:
    remainder = decimal_number % hex_base
    hexadecimal_number = hex_digits[remainder] + hexadecimal_number
    decimal_number //= hex_base

print(hexadecimal_number)