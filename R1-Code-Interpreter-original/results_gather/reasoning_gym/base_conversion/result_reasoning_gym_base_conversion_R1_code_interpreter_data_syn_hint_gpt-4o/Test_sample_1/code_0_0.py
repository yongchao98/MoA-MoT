# Step 1: Convert base-4 number to decimal
base_4_number = "21020"
decimal_number = 0
base_4_length = len(base_4_number)

for i in range(base_4_length):
    digit = int(base_4_number[base_4_length - i - 1])
    decimal_number += digit * (4 ** i)

# Step 2: Convert decimal number to base-12
base_12_number = ""
while decimal_number > 0:
    remainder = decimal_number % 12
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_12_number = chr(remainder - 10 + ord('a')) + base_12_number
    else:
        base_12_number = str(remainder) + base_12_number
    decimal_number //= 12

print(base_12_number)