# Step 1: Convert base-5 to decimal
base_5_number = "12141"
decimal_number = 0
base_5_length = len(base_5_number)

for i in range(base_5_length):
    digit = int(base_5_number[i])
    power = base_5_length - i - 1
    decimal_number += digit * (5 ** power)

# Step 2: Convert decimal to base-15
base_15_number = ""
while decimal_number > 0:
    remainder = decimal_number % 15
    if remainder >= 10:
        # Convert remainder to corresponding letter
        base_15_number = chr(remainder - 10 + ord('a')) + base_15_number
    else:
        base_15_number = str(remainder) + base_15_number
    decimal_number //= 15

print(base_15_number)