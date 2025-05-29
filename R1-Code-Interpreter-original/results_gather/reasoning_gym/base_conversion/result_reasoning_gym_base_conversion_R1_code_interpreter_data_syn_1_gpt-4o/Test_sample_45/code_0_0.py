# Step 1: Convert base-6 number 3442 to decimal
base_6_number = "3442"
base_6 = 6
decimal_number = 0

# Calculate the decimal equivalent
for i, digit in enumerate(reversed(base_6_number)):
    decimal_number += int(digit) * (base_6 ** i)

# Step 2: Convert the decimal number to base-7
base_7 = 7
base_7_number = ""

# Calculate the base-7 equivalent
while decimal_number > 0:
    remainder = decimal_number % base_7
    base_7_number = str(remainder) + base_7_number
    decimal_number //= base_7

print(base_7_number)