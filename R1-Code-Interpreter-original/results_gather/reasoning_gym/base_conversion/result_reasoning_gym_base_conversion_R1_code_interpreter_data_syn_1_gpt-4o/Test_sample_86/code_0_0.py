# Step 1: Convert base-12 number 508 to base-10
base_12_number = "508"
base_12 = 12

# Convert base-12 to base-10
base_10_number = 0
for i, digit in enumerate(reversed(base_12_number)):
    base_10_number += int(digit) * (base_12 ** i)

# Step 2: Convert base-10 number to base-4
base_4 = 4
base_4_number = ""

while base_10_number > 0:
    remainder = base_10_number % base_4
    base_4_number = str(remainder) + base_4_number
    base_10_number //= base_4

print(base_4_number)