# Step 1: Convert base-11 number 748 to base-10
base_11_number = "748"
base_11 = 11

# Convert to base-10
base_10_number = 0
for index, digit in enumerate(reversed(base_11_number)):
    base_10_number += int(digit) * (base_11 ** index)

# Step 2: Convert base-10 number to base-7
base_7 = 7
base_7_number = ""

while base_10_number > 0:
    remainder = base_10_number % base_7
    base_7_number = str(remainder) + base_7_number
    base_10_number //= base_7

print(base_7_number)