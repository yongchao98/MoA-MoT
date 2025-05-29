# Step 1: Convert base-9 number 174 to decimal
base_9_number = "174"
decimal_number = 0
base_9 = 9

# Calculate the decimal equivalent
for i, digit in enumerate(reversed(base_9_number)):
    decimal_number += int(digit) * (base_9 ** i)

# Step 2: Convert the decimal number to base-7
base_7 = 7
base_7_number = ""

while decimal_number > 0:
    remainder = decimal_number % base_7
    base_7_number = str(remainder) + base_7_number
    decimal_number //= base_7

print(base_7_number)