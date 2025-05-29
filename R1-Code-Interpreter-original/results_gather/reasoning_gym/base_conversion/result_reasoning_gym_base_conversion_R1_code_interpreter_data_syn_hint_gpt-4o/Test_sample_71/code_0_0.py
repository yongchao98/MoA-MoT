# Convert base-9 number 1002 to decimal
base_9_number = "1002"
decimal_number = 0
base_9_length = len(base_9_number)

for i in range(base_9_length):
    digit = int(base_9_number[base_9_length - i - 1])
    decimal_number += digit * (9 ** i)

# Convert decimal number to hexadecimal
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