# Step 1: Convert base-7 to decimal
base_7_number = "1431"
decimal_number = 0
base_7_length = len(base_7_number)

for i in range(base_7_length):
    digit = int(base_7_number[i])
    power = base_7_length - i - 1
    decimal_number += digit * (7 ** power)

# Step 2: Convert decimal to base-14
base_14_number = ""
while decimal_number > 0:
    remainder = decimal_number % 14
    if remainder >= 10:
        # Convert remainder to corresponding letter for base > 10
        base_14_number = chr(remainder - 10 + ord('a')) + base_14_number
    else:
        base_14_number = str(remainder) + base_14_number
    decimal_number //= 14

print(base_14_number)