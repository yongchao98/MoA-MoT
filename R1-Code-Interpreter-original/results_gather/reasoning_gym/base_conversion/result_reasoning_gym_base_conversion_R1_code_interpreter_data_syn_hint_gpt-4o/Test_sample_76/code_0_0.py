# Convert base-13 number 'a3' to decimal
base_13_number = 'a3'
base_13 = 13

# Convert to decimal
decimal_value = 0
for index, digit in enumerate(reversed(base_13_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_13 ** index)

# Convert decimal to base-15
base_15 = 15
base_15_number = ''
while decimal_value > 0:
    remainder = decimal_value % base_15
    if remainder < 10:
        base_15_number = str(remainder) + base_15_number
    else:
        base_15_number = chr(remainder - 10 + ord('a')) + base_15_number
    decimal_value //= base_15

print(base_15_number)