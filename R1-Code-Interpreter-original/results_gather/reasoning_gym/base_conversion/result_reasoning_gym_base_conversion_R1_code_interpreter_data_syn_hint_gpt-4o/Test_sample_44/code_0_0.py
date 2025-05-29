# Step 1: Convert base-9 number 461 to decimal
base_9_number = '461'
decimal_number = 0
base_9 = 9

# Calculate the decimal equivalent
for index, digit in enumerate(reversed(base_9_number)):
    decimal_number += int(digit) * (base_9 ** index)

# Step 2: Convert the decimal number to base-8
base_8 = 8
base_8_number = ''

while decimal_number > 0:
    remainder = decimal_number % base_8
    base_8_number = str(remainder) + base_8_number
    decimal_number //= base_8

print(base_8_number)