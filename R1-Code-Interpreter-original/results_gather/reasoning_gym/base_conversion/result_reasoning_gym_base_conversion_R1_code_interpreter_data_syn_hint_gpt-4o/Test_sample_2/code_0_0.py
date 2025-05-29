# Step 1: Convert base-7 number 1402 to decimal
base_7_number = "1402"
decimal_number = 0
base_7_length = len(base_7_number)

for i in range(base_7_length):
    digit = int(base_7_number[i])
    power = base_7_length - i - 1
    decimal_number += digit * (7 ** power)

# Step 2: Convert the decimal number to base-11
base_11_number = ""
while decimal_number > 0:
    remainder = decimal_number % 11
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_11_number = chr(remainder - 10 + ord('a')) + base_11_number
    else:
        base_11_number = str(remainder) + base_11_number
    decimal_number //= 11

print(base_11_number)