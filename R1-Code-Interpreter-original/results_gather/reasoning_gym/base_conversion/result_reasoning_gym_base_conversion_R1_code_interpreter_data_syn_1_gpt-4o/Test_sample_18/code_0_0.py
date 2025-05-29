# Step 1: Convert from base-15 to decimal
base_15_number = "23e"
base_15 = 15

# Convert each digit to its decimal equivalent
decimal_value = 0
for i, digit in enumerate(reversed(base_15_number)):
    if digit.isdigit():
        value = int(digit)
    else:
        value = ord(digit) - ord('a') + 10
    decimal_value += value * (base_15 ** i)

# Step 2: Convert from decimal to base-6
base_6 = 6
base_6_number = ""

while decimal_value > 0:
    remainder = decimal_value % base_6
    base_6_number = str(remainder) + base_6_number
    decimal_value //= base_6

print(base_6_number)