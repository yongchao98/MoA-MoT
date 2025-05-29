# Step 1: Convert base-3 to base-10
base_3_number = "2211"
base_10_number = 0
for i, digit in enumerate(reversed(base_3_number)):
    base_10_number += int(digit) * (3 ** i)

# Step 2: Convert base-10 to base-14
base_14_number = ""
while base_10_number > 0:
    remainder = base_10_number % 14
    if remainder >= 10:
        # Convert remainder to corresponding letter
        base_14_number = chr(remainder - 10 + ord('a')) + base_14_number
    else:
        base_14_number = str(remainder) + base_14_number
    base_10_number //= 14

print(base_14_number)