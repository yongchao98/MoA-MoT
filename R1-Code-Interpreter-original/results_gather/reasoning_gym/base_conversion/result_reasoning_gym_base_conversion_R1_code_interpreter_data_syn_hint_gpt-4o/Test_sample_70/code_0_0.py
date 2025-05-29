# Step 1: Convert base-12 number 107 to base-10
base_12_number = "107"
base_12 = 12

# Convert to base-10
base_10_number = 0
for i, digit in enumerate(reversed(base_12_number)):
    base_10_number += int(digit) * (base_12 ** i)

# Step 2: Convert base-10 number to base-11
base_11 = 11
base_11_number = ""

while base_10_number > 0:
    remainder = base_10_number % base_11
    base_11_number = str(remainder) + base_11_number
    base_10_number //= base_11

print(base_11_number)