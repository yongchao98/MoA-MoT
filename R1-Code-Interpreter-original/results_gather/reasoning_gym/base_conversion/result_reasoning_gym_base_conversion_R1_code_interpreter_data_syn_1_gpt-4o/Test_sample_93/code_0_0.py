# Step 1: Convert base-5 to decimal
base_5_number = "12141"
decimal_number = 0
base_5_number_reversed = base_5_number[::-1]

for i, digit in enumerate(base_5_number_reversed):
    decimal_number += int(digit) * (5 ** i)

# Step 2: Convert decimal to base-15
base_15_number = ""
while decimal_number > 0:
    remainder = decimal_number % 15
    if remainder >= 10:
        # Convert remainder to corresponding letter
        base_15_number = chr(remainder - 10 + ord('a')) + base_15_number
    else:
        base_15_number = str(remainder) + base_15_number
    decimal_number //= 15

print(base_15_number)