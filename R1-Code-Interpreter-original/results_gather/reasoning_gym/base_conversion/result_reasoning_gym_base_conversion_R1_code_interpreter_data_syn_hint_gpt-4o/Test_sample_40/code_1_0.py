# Step 1: Convert base-14 number "55" to decimal
base_14_number = "55"
base_14 = 14
decimal_value = 0

# Calculate the decimal value
for i, digit in enumerate(reversed(base_14_number)):
    decimal_value += int(digit) * (base_14 ** i)

# Step 2: Convert the decimal value to base-11
base_11 = 11
base_11_number = ""

while decimal_value > 0:
    remainder = decimal_value % base_11
    base_11_number = str(remainder) + base_11_number
    decimal_value //= base_11

print(base_11_number)