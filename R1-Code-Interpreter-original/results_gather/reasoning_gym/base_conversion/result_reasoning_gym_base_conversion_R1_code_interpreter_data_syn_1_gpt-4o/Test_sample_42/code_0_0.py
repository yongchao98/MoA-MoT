# Base-11 number
base_11_number = "730"

# Convert each digit to base-10
base_10_value = 0
base = 11

# Calculate the base-10 equivalent
for i, digit in enumerate(reversed(base_11_number)):
    base_10_value += int(digit) * (base ** i)

print(base_10_value)