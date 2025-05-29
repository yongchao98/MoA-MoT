# Step 1: Convert base-15 number "23d" to base-10
base_15_number = "23d"
base_15_digits = "0123456789abcde"
base_10_value = 0

# Calculate base-10 equivalent
for i, digit in enumerate(reversed(base_15_number)):
    base_10_value += base_15_digits.index(digit) * (15 ** i)

# Step 2: Convert base-10 number to base-9
base_9_value = ""
while base_10_value > 0:
    remainder = base_10_value % 9
    base_9_value = str(remainder) + base_9_value
    base_10_value //= 9

print(base_9_value)