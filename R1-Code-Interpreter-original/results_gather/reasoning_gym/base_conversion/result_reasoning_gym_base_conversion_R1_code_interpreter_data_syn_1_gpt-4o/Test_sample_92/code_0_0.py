# Base-6 number
base_6_number = "4315"

# Convert base-6 to base-10
base_10_number = sum(int(digit) * (6 ** idx) for idx, digit in enumerate(reversed(base_6_number)))

print(base_10_number)