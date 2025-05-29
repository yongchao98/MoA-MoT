# Step 1: Convert base-5 to base-10
base_5_number = "4334"
base_10_number = 0
base_5_number_reversed = base_5_number[::-1]

for i, digit in enumerate(base_5_number_reversed):
    base_10_number += int(digit) * (5 ** i)

# Step 2: Convert base-10 to base-13
base_13_number = ""
current_number = base_10_number

while current_number > 0:
    remainder = current_number % 13
    if remainder >= 10:
        # Convert remainder to corresponding letter
        base_13_number = chr(remainder - 10 + ord('a')) + base_13_number
    else:
        base_13_number = str(remainder) + base_13_number
    current_number //= 13

print(base_13_number)