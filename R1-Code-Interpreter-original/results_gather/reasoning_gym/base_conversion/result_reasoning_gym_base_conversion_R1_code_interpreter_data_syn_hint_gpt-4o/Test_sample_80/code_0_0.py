# Convert base-15 number "22c" to decimal
base_15_number = "22c"
base_15 = 15

# Convert to decimal
decimal_value = 0
for i, digit in enumerate(reversed(base_15_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_15 ** i)

# Convert decimal to base-12
base_12 = 12
base_12_number = ""
while decimal_value > 0:
    remainder = decimal_value % base_12
    if remainder >= 10:
        base_12_number = chr(remainder - 10 + ord('a')) + base_12_number
    else:
        base_12_number = str(remainder) + base_12_number
    decimal_value //= base_12

print(base_12_number)