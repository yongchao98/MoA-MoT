# Step 1: Convert base-15 number 424 to base-10
base_15_number = "424"
base_15 = 15

# Convert to base-10
base_10_number = 0
for i, digit in enumerate(reversed(base_15_number)):
    base_10_number += int(digit) * (base_15 ** i)

# Step 2: Convert base-10 number to base-9
base_9 = 9
base_9_number = ""

while base_10_number > 0:
    remainder = base_10_number % base_9
    base_9_number = str(remainder) + base_9_number
    base_10_number //= base_9

print(base_9_number)