# Step 1: Convert base-15 to base-10
base_15_number = '442'
base_15 = 15

# Convert to decimal
decimal_number = 0
for i, digit in enumerate(reversed(base_15_number)):
    decimal_number += int(digit, base_15) * (base_15 ** i)

# Step 2: Convert base-10 to base-5
base_5 = 5
base_5_number = ''

while decimal_number > 0:
    remainder = decimal_number % base_5
    base_5_number = str(remainder) + base_5_number
    decimal_number //= base_5

print(base_5_number)