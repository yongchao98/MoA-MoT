# Step 1: Convert base-11 number 204 to base-10
base_11_number = "204"
base_11 = 11

# Convert to base-10
base_10_number = 0
for i, digit in enumerate(reversed(base_11_number)):
    base_10_number += int(digit) * (base_11 ** i)

# Step 2: Convert base-10 number to base-3
base_3 = 3
base_3_number = ""

while base_10_number > 0:
    remainder = base_10_number % base_3
    base_3_number = str(remainder) + base_3_number
    base_10_number //= base_3

print(base_3_number)