# Step 1: Convert base-6 number 103 to decimal
base_6_number = "103"
decimal_number = 0
base_6_number_reversed = base_6_number[::-1]

for i, digit in enumerate(base_6_number_reversed):
    decimal_number += int(digit) * (6 ** i)

# Step 2: Convert the decimal number to base-13
base_13_number = ""
current_number = decimal_number

while current_number > 0:
    remainder = current_number % 13
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_13_number = chr(remainder - 10 + ord('a')) + base_13_number
    else:
        base_13_number = str(remainder) + base_13_number
    current_number //= 13

print(base_13_number)