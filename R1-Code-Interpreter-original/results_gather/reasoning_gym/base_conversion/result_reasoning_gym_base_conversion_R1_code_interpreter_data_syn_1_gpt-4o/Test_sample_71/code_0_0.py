# Step 1: Convert base-9 number to decimal
base_9_number = "1002"
decimal_number = 0

# Convert base-9 to decimal
for index, digit in enumerate(reversed(base_9_number)):
    decimal_number += int(digit) * (9 ** index)

# Step 2: Convert decimal number to hexadecimal
hexadecimal_number = ""
hex_digits = "0123456789abcdef"

while decimal_number > 0:
    remainder = decimal_number % 16
    hexadecimal_number = hex_digits[remainder] + hexadecimal_number
    decimal_number //= 16

# If the number is 0, we should return "0"
if not hexadecimal_number:
    hexadecimal_number = "0"

print(hexadecimal_number)