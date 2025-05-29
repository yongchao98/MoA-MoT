# Step 1: Convert base-11 number 204 to decimal
base_11_number = "204"
decimal_number = 0
base_11 = 11

# Calculate the decimal equivalent
for i, digit in enumerate(reversed(base_11_number)):
    decimal_number += int(digit) * (base_11 ** i)

# Step 2: Convert the decimal number to base-3
base_3 = 3
base_3_number = ""

while decimal_number > 0:
    remainder = decimal_number % base_3
    base_3_number = str(remainder) + base_3_number
    decimal_number //= base_3

print(base_3_number)