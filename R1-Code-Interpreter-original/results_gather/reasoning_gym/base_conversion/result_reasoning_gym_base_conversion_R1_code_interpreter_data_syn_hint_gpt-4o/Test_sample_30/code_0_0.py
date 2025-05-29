# Step 1: Convert base-11 number 794 to decimal
base_11_number = "794"
base_11 = 11

# Convert to decimal
decimal_number = 0
for index, digit in enumerate(reversed(base_11_number)):
    decimal_number += int(digit) * (base_11 ** index)

# Step 2: Convert decimal number to base-4
base_4 = 4
base_4_number = ""

while decimal_number > 0:
    remainder = decimal_number % base_4
    base_4_number = str(remainder) + base_4_number
    decimal_number //= base_4

print(base_4_number)