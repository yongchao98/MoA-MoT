# Step 1: Convert base-5 to decimal
base_5_number = "4102"
decimal_number = 0
base = 5

# Convert base-5 to decimal
for i, digit in enumerate(reversed(base_5_number)):
    decimal_number += int(digit) * (base ** i)

# Step 2: Convert decimal to hexadecimal
hexadecimal_number = ""
hex_digits = "0123456789abcdef"

while decimal_number > 0:
    remainder = decimal_number % 16
    hexadecimal_number = hex_digits[remainder] + hexadecimal_number
    decimal_number //= 16

# If the number is 0, the hexadecimal representation is also 0
if not hexadecimal_number:
    hexadecimal_number = "0"

print(hexadecimal_number)