# Step 1: Convert base-11 number 301 to base-10
base_11_number = "301"
base_11 = 11

# Convert to base-10
base_10_number = 0
for i, digit in enumerate(reversed(base_11_number)):
    base_10_number += int(digit) * (base_11 ** i)

# Step 2: Convert base-10 number to base-15
base_15 = 15
base_15_number = ""

while base_10_number > 0:
    remainder = base_10_number % base_15
    if remainder >= 10:
        # Convert remainder to corresponding letter for bases > 10
        base_15_number = chr(remainder - 10 + ord('a')) + base_15_number
    else:
        base_15_number = str(remainder) + base_15_number
    base_10_number //= base_15

print(base_15_number)