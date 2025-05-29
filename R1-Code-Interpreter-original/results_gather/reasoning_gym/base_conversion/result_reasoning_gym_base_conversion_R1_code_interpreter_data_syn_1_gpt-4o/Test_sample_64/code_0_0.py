# Step 1: Convert base-7 number to decimal
base_7_number = "1051"
decimal_number = 0
base = 7

# Convert base-7 to decimal
for index, digit in enumerate(reversed(base_7_number)):
    decimal_number += int(digit) * (base ** index)

# Step 2: Convert decimal number to hexadecimal
hexadecimal_number = ""
hex_base = 16
hex_digits = "0123456789abcdef"

while decimal_number > 0:
    remainder = decimal_number % hex_base
    hexadecimal_number = hex_digits[remainder] + hexadecimal_number
    decimal_number //= hex_base

# If the number is 0, the hexadecimal representation should be "0"
if not hexadecimal_number:
    hexadecimal_number = "0"

print(hexadecimal_number)