# Step 1: Convert base-5 number 4042 to base-10
base_5_number = "4042"
base_5 = 5

# Convert base-5 to base-10
base_10_number = 0
for index, digit in enumerate(reversed(base_5_number)):
    base_10_number += int(digit) * (base_5 ** index)

# Step 2: Convert base-10 number to base-11
base_11 = 11
base_11_number = ""

while base_10_number > 0:
    remainder = base_10_number % base_11
    if remainder >= 10:
        # Convert remainder to corresponding letter for bases > 10
        base_11_number = chr(remainder - 10 + ord('a')) + base_11_number
    else:
        base_11_number = str(remainder) + base_11_number
    base_10_number //= base_11

print(base_11_number)