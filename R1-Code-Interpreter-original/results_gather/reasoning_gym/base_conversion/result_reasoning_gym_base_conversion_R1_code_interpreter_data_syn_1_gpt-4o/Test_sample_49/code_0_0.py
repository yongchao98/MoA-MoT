# Step 1: Convert from base-15 to base-10
base_15_number = "248"
base_15 = 15

# Convert base-15 to decimal
decimal_number = 0
for index, digit in enumerate(reversed(base_15_number)):
    decimal_number += int(digit) * (base_15 ** index)

# Step 2: Convert from decimal to base-9
base_9 = 9
base_9_number = ""

while decimal_number > 0:
    remainder = decimal_number % base_9
    base_9_number = str(remainder) + base_9_number
    decimal_number //= base_9

print(base_9_number)