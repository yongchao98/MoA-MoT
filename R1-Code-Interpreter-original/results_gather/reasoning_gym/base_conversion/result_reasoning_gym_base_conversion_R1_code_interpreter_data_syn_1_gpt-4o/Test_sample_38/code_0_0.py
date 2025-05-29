# Convert base-3 number to base-10
base_3_number = '1021111'
base_10_number = 0

# Calculate base-10 equivalent
for index, digit in enumerate(reversed(base_3_number)):
    base_10_number += int(digit) * (3 ** index)

# Convert base-10 number to base-6
base_6_number = ''
while base_10_number > 0:
    remainder = base_10_number % 6
    base_6_number = str(remainder) + base_6_number
    base_10_number //= 6

print(base_6_number)