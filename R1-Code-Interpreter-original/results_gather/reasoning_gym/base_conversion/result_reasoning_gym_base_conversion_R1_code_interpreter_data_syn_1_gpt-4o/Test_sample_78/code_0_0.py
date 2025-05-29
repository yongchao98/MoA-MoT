# Step 1: Convert base-8 number to base-10
base_8_number = '1674'
base_10_number = 0

# Convert base-8 to base-10
for i, digit in enumerate(reversed(base_8_number)):
    base_10_number += int(digit) * (8 ** i)

# Step 2: Convert base-10 number to base-14
base_14_number = ''
digits = '0123456789abcdef'

while base_10_number > 0:
    remainder = base_10_number % 14
    base_14_number = digits[remainder] + base_14_number
    base_10_number //= 14

print(base_14_number)