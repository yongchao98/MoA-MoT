# Step 1: Convert base-11 number "20a" to decimal
base_11_number = "20a"
base_11 = 11

# Convert each digit to its decimal equivalent
decimal_value = 0
for i, digit in enumerate(reversed(base_11_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_11 ** i)

# Step 2: Convert the decimal number to base-14
base_14 = 14
base_14_number = ""

while decimal_value > 0:
    remainder = decimal_value % base_14
    if remainder < 10:
        base_14_number = str(remainder) + base_14_number
    else:
        base_14_number = chr(remainder - 10 + ord('a')) + base_14_number
    decimal_value //= base_14

print(base_14_number)