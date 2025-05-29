# Convert base-4 number 3311 to base-10
base_4_number = '3311'
base_10_number = 0
base_4 = 4

# Convert base-4 to base-10
for index, digit in enumerate(reversed(base_4_number)):
    base_10_number += int(digit) * (base_4 ** index)

# Convert base-10 to base-5
base_5_number = ''
base_5 = 5

while base_10_number > 0:
    remainder = base_10_number % base_5
    base_5_number = str(remainder) + base_5_number
    base_10_number //= base_5

print(base_5_number)