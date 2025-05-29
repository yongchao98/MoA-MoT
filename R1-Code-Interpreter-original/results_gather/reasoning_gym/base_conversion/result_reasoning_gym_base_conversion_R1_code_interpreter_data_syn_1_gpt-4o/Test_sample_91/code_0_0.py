# Step 1: Convert base-6 number 125 to base-10
base_6_number = "125"
base_10_number = 0
base_6_number_reversed = base_6_number[::-1]

for i, digit in enumerate(base_6_number_reversed):
    base_10_number += int(digit) * (6 ** i)

# Step 2: Convert base-10 number to base-3
base_3_number = ""
while base_10_number > 0:
    remainder = base_10_number % 3
    base_3_number = str(remainder) + base_3_number
    base_10_number //= 3

print(base_3_number)