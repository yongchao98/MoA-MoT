# Step 1: Convert base-15 number '2a0' to decimal
base_15_number = '2a0'
base_15 = 15

# Convert each digit to decimal
decimal_value = 0
for i, digit in enumerate(reversed(base_15_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_15 ** i)

# Step 2: Convert the decimal number to base-14
base_14 = 14
base_14_number = ''

while decimal_value > 0:
    remainder = decimal_value % base_14
    if remainder < 10:
        base_14_number = str(remainder) + base_14_number
    else:
        base_14_number = chr(remainder - 10 + ord('a')) + base_14_number
    decimal_value //= base_14

print(base_14_number)