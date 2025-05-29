# Step 1: Convert base-6 to base-10
base_6_number = "4220"
base_10_number = 0
base_6_number_reversed = base_6_number[::-1]

for i, digit in enumerate(base_6_number_reversed):
    base_10_number += int(digit) * (6 ** i)

# Step 2: Convert base-10 to base-4
base_4_number = ""
while base_10_number > 0:
    remainder = base_10_number % 4
    base_4_number = str(remainder) + base_4_number
    base_10_number //= 4

print(base_4_number)