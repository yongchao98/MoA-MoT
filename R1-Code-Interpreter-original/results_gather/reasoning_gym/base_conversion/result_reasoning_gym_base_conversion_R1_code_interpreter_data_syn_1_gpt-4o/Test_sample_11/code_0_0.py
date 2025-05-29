# Step 1: Convert base-6 number 103 to base-10
base_6_number = "103"
base_10_number = 0
base_6_length = len(base_6_number)

for i in range(base_6_length):
    digit = int(base_6_number[i])
    power = base_6_length - i - 1
    base_10_number += digit * (6 ** power)

# Step 2: Convert base-10 number to base-13
base_13_number = ""
base = 13
while base_10_number > 0:
    remainder = base_10_number % base
    if remainder >= 10:
        # Convert remainder to corresponding letter
        remainder = chr(remainder - 10 + ord('a'))
    base_13_number = str(remainder) + base_13_number
    base_10_number //= base

print(base_13_number)