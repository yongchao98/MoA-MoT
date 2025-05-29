# Convert base-7 number 1431 to base-10
base_7_number = "1431"
base_10_number = 0
base_7 = 7

# Convert base-7 to base-10
for index, digit in enumerate(reversed(base_7_number)):
    base_10_number += int(digit) * (base_7 ** index)

# Convert base-10 to base-4
base_4_number = ""
base_4 = 4

while base_10_number > 0:
    remainder = base_10_number % base_4
    base_4_number = str(remainder) + base_4_number
    base_10_number //= base_4

print(base_4_number)