# Step 1: Convert base-8 number to decimal
base_8_number = "1674"
decimal_number = 0
for i, digit in enumerate(reversed(base_8_number)):
    decimal_number += int(digit) * (8 ** i)

# Step 2: Convert decimal number to base-14
base_14_number = ""
while decimal_number > 0:
    remainder = decimal_number % 14
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_14_number = chr(remainder - 10 + ord('a')) + base_14_number
    else:
        base_14_number = str(remainder) + base_14_number
    decimal_number //= 14

print(base_14_number)