# Step 1: Convert base-13 number "130" to decimal
base_13_number = "130"
base_13 = 13

# Calculate the decimal equivalent
decimal_value = 0
for i, digit in enumerate(reversed(base_13_number)):
    decimal_value += int(digit) * (base_13 ** i)

# Step 2: Convert the decimal number to base-3
base_3 = 3
base_3_number = ""

while decimal_value > 0:
    remainder = decimal_value % base_3
    base_3_number = str(remainder) + base_3_number
    decimal_value //= base_3

print(base_3_number)