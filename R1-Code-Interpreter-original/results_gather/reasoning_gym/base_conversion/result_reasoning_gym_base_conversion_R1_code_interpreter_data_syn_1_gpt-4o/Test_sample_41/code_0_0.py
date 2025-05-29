# Base-3 number
base_3_number = "101121"

# Convert base-3 to base-10
base_10_number = 0
for i, digit in enumerate(reversed(base_3_number)):
    base_10_number += int(digit) * (3 ** i)

print(base_10_number)