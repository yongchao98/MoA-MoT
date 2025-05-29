# Step 1: Convert base-6 number to decimal
base_6_number = '3245'
decimal_number = 0
base = 6

# Convert base-6 to decimal
for index, digit in enumerate(reversed(base_6_number)):
    decimal_number += int(digit) * (base ** index)

# Step 2: Convert decimal to hexadecimal
hexadecimal_number = ''
hex_base = 16
hex_digits = '0123456789abcdef'

while decimal_number > 0:
    remainder = decimal_number % hex_base
    hexadecimal_number = hex_digits[remainder] + hexadecimal_number
    decimal_number //= hex_base

print(hexadecimal_number)