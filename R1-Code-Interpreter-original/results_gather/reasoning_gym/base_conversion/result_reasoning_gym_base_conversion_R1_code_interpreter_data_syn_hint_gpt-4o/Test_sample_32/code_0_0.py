# Step 1: Convert base-3 number to decimal
base_3_number = "200001"
decimal_number = 0

# Convert base-3 to decimal
for index, digit in enumerate(reversed(base_3_number)):
    decimal_number += int(digit) * (3 ** index)

# Step 2: Convert decimal to base-12
base_12_number = ""
while decimal_number > 0:
    remainder = decimal_number % 12
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_12_number = chr(remainder - 10 + ord('a')) + base_12_number
    else:
        base_12_number = str(remainder) + base_12_number
    decimal_number //= 12

print(base_12_number)