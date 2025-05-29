# Step 1: Convert base-9 number 461 to decimal
base_9_number = "461"
decimal_number = 0
base_9 = 9

# Convert base-9 to decimal
for i, digit in enumerate(reversed(base_9_number)):
    decimal_number += int(digit) * (base_9 ** i)

# Step 2: Convert decimal number to base-8
base_8 = 8
base_8_number = ""

# Convert decimal to base-8
while decimal_number > 0:
    remainder = decimal_number % base_8
    base_8_number = str(remainder) + base_8_number
    decimal_number //= base_8

print(base_8_number)