# Step 1: Convert base-14 to decimal
base_14_number = "5d"
base_14 = 14

# Convert each digit to decimal
decimal_value = 0
for i, digit in enumerate(reversed(base_14_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_14 ** i)

# Step 2: Convert decimal to base-11
base_11 = 11
base_11_number = ""

while decimal_value > 0:
    remainder = decimal_value % base_11
    if remainder >= 10:
        base_11_number = chr(remainder - 10 + ord('a')) + base_11_number
    else:
        base_11_number = str(remainder) + base_11_number
    decimal_value //= base_11

print(base_11_number)