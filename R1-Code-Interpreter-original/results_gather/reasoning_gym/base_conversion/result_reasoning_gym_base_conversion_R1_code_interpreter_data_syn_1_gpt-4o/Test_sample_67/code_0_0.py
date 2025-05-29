# Step 1: Convert base-8 number 635 to base-10
base_8_number = '635'
base_10_number = 0
base_8 = 8

# Convert base-8 to base-10
for index, digit in enumerate(reversed(base_8_number)):
    base_10_number += int(digit) * (base_8 ** index)

# Step 2: Convert base-10 number to base-5
base_5_number = ''
base_5 = 5

# Convert base-10 to base-5
while base_10_number > 0:
    remainder = base_10_number % base_5
    base_5_number = str(remainder) + base_5_number
    base_10_number //= base_5

print(base_5_number)