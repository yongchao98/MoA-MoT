# Step 1: Convert base-5 number 332 to decimal
base_5_number = '332'
decimal_number = 0
base_5 = 5

# Calculate the decimal equivalent
for index, digit in enumerate(reversed(base_5_number)):
    decimal_number += int(digit) * (base_5 ** index)

# Step 2: Convert the decimal number to base-4
base_4 = 4
base_4_number = ''

# Calculate the base-4 equivalent
while decimal_number > 0:
    remainder = decimal_number % base_4
    base_4_number = str(remainder) + base_4_number
    decimal_number //= base_4

print(base_4_number)