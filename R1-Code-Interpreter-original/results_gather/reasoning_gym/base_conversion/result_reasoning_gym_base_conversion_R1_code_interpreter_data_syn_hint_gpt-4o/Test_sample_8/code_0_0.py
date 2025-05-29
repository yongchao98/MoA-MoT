# Step 1: Convert base-11 number to decimal
base_11_number = "397"
decimal_number = 0
base = 11

# Calculate the decimal equivalent
for i, digit in enumerate(reversed(base_11_number)):
    decimal_number += int(digit) * (base ** i)

# Step 2: Convert decimal number to hexadecimal
hexadecimal_number = ""
hex_digits = "0123456789abcdef"

while decimal_number > 0:
    remainder = decimal_number % 16
    hexadecimal_number = hex_digits[remainder] + hexadecimal_number
    decimal_number //= 16

print(hexadecimal_number)