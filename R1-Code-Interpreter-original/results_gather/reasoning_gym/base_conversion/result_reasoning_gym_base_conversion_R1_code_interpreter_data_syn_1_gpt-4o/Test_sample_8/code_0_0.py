# Step 1: Convert base-11 number to decimal
base_11_number = "397"
decimal_number = 0
base_11_length = len(base_11_number)

for i, digit in enumerate(reversed(base_11_number)):
    decimal_number += int(digit) * (11 ** i)

# Step 2: Convert decimal number to hexadecimal
hexadecimal_number = ""
while decimal_number > 0:
    remainder = decimal_number % 16
    if remainder < 10:
        hexadecimal_number = str(remainder) + hexadecimal_number
    else:
        # Convert remainder to corresponding hexadecimal character
        hexadecimal_number = chr(ord('a') + remainder - 10) + hexadecimal_number
    decimal_number //= 16

print(hexadecimal_number)