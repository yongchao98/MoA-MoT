# Step 1: Convert base-4 to base-10
base_4_number = "21020"
base_10_number = 0

# Calculate the base-10 equivalent
for index, digit in enumerate(reversed(base_4_number)):
    base_10_number += int(digit) * (4 ** index)

# Step 2: Convert base-10 to base-12
base_12_number = ""
while base_10_number > 0:
    remainder = base_10_number % 12
    if remainder >= 10:
        # Convert remainder to a letter if it's 10 or 11
        base_12_number = chr(remainder - 10 + ord('a')) + base_12_number
    else:
        base_12_number = str(remainder) + base_12_number
    base_10_number //= 12

print(base_12_number)