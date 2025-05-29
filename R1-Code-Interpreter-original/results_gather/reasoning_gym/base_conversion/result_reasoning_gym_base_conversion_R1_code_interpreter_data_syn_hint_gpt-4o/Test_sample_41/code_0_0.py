# Base-3 number
base_3_number = "101121"

# Convert base-3 to base-10
base_10_value = 0
for index, digit in enumerate(reversed(base_3_number)):
    base_10_value += int(digit) * (3 ** index)

print(base_10_value)