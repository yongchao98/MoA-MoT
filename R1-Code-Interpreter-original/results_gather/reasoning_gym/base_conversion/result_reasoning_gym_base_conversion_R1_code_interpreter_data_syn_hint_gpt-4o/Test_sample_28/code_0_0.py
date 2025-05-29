# Convert base-6 number 1243 to base-10
base_6_number = "1243"
base_6 = 6

# Convert to decimal
decimal_number = 0
for index, digit in enumerate(reversed(base_6_number)):
    decimal_number += int(digit) * (base_6 ** index)

# Convert decimal to base-11
base_11 = 11
base_11_number = ""
while decimal_number > 0:
    remainder = decimal_number % base_11
    if remainder >= 10:
        # Convert remainder to corresponding letter
        base_11_number = chr(remainder - 10 + ord('a')) + base_11_number
    else:
        base_11_number = str(remainder) + base_11_number
    decimal_number //= base_11

print(base_11_number)