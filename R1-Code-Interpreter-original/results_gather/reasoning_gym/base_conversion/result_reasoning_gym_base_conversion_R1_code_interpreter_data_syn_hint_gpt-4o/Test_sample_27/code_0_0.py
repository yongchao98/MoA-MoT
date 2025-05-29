# Step 1: Convert base-3 number to decimal
base_3_number = '22111'
decimal_number = 0
base = 3

# Convert base-3 to decimal
for i, digit in enumerate(reversed(base_3_number)):
    decimal_number += int(digit) * (base ** i)

# Step 2: Convert decimal to hexadecimal
hexadecimal_number = ''
hex_digits = '0123456789abcdef'

while decimal_number > 0:
    remainder = decimal_number % 16
    hexadecimal_number = hex_digits[remainder] + hexadecimal_number
    decimal_number //= 16

print(hexadecimal_number)