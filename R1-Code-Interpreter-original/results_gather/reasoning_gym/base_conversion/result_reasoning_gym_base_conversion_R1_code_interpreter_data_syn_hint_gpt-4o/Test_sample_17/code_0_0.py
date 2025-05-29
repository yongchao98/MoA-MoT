# Step 1: Convert base-7 to base-10
base_7_number = "125"
base_10_number = 0
for i, digit in enumerate(reversed(base_7_number)):
    base_10_number += int(digit) * (7 ** i)

# Step 2: Convert base-10 to base-8
base_8_number = ""
while base_10_number > 0:
    remainder = base_10_number % 8
    base_8_number = str(remainder) + base_8_number
    base_10_number //= 8

print(base_8_number)