# Step 1: Convert base-3 number to decimal
base_3_number = '2211'
decimal_number = 0
base_3_length = len(base_3_number)

for i in range(base_3_length):
    digit = int(base_3_number[base_3_length - i - 1])
    decimal_number += digit * (3 ** i)

# Step 2: Convert decimal number to base-14
base_14_number = ''
while decimal_number > 0:
    remainder = decimal_number % 14
    if remainder >= 10:
        # Convert remainder to corresponding letter
        base_14_number = chr(remainder - 10 + ord('a')) + base_14_number
    else:
        base_14_number = str(remainder) + base_14_number
    decimal_number //= 14

print(base_14_number)