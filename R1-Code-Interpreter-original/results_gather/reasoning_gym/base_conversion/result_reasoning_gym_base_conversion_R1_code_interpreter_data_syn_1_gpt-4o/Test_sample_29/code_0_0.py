# Base-4 number
base_4_number = "13322"

# Convert base-4 to base-10
base_10_number = sum(int(digit) * (4 ** idx) for idx, digit in enumerate(reversed(base_4_number)))

print(base_10_number)