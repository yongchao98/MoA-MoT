# Convert base-6 number 4125 to base-10
base_6_number = "4125"
base_10_number = 0
base_6_number_reversed = base_6_number[::-1]

for i, digit in enumerate(base_6_number_reversed):
    base_10_number += int(digit) * (6 ** i)

# Convert base-10 number to base-5
base_5_number = ""
while base_10_number > 0:
    remainder = base_10_number % 5
    base_5_number = str(remainder) + base_5_number
    base_10_number //= 5

print(base_5_number)