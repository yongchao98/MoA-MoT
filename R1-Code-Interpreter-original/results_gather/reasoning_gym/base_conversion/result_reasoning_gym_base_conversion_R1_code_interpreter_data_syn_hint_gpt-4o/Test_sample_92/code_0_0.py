# Base-6 number
base_6_number = "4315"

# Convert to base-10
base_10_number = 0
for i, digit in enumerate(reversed(base_6_number)):
    base_10_number += int(digit) * (6 ** i)

print(base_10_number)