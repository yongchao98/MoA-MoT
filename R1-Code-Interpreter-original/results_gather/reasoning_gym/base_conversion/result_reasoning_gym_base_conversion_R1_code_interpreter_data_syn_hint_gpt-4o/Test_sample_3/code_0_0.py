# Convert base-7 number 1431 to base-10
base_7_number = "1431"
base_10_number = 0
base_7_number_reversed = base_7_number[::-1]

for i, digit in enumerate(base_7_number_reversed):
    base_10_number += int(digit) * (7 ** i)

# Convert base-10 number to base-14
base_14_number = ""
current_number = base_10_number

while current_number > 0:
    remainder = current_number % 14
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_14_number = chr(remainder - 10 + ord('a')) + base_14_number
    else:
        base_14_number = str(remainder) + base_14_number
    current_number //= 14

print(base_14_number)