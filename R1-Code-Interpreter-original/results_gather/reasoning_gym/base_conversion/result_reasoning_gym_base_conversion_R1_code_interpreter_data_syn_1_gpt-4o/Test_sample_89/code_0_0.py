# Step 1: Convert base-5 number 324 to base-10
base_5_number = '324'
base_5 = 5
decimal_number = 0

# Calculate the decimal equivalent
for i, digit in enumerate(reversed(base_5_number)):
    decimal_number += int(digit) * (base_5 ** i)

# Step 2: Convert the decimal number to base-11
base_11 = 11
base_11_number = ''

# Convert to base-11
while decimal_number > 0:
    remainder = decimal_number % base_11
    if remainder >= 10:
        # Convert remainder to a letter if it's 10 or more
        base_11_number = chr(remainder - 10 + ord('a')) + base_11_number
    else:
        base_11_number = str(remainder) + base_11_number
    decimal_number //= base_11

print(base_11_number)