# Step 1: Convert base-12 number '27' to base-10
base_12_number = '27'
base_10_number = 0
base_12 = 12

# Convert each digit to base-10
for i, digit in enumerate(reversed(base_12_number)):
    base_10_number += int(digit) * (base_12 ** i)

# Step 2: Convert base-10 number to base-15
base_15_number = ''
base_15 = 15

while base_10_number > 0:
    remainder = base_10_number % base_15
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_15_number = chr(remainder - 10 + ord('a')) + base_15_number
    else:
        base_15_number = str(remainder) + base_15_number
    base_10_number //= base_15

print(base_15_number)