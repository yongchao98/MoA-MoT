# Step 1: Convert base-3 to base-10
base_3_number = '1020200'
base_10_number = 0
for i, digit in enumerate(reversed(base_3_number)):
    base_10_number += int(digit) * (3 ** i)

# Step 2: Convert base-10 to base-8
base_8_number = ''
while base_10_number > 0:
    remainder = base_10_number % 8
    base_8_number = str(remainder) + base_8_number
    base_10_number //= 8

print(base_8_number)