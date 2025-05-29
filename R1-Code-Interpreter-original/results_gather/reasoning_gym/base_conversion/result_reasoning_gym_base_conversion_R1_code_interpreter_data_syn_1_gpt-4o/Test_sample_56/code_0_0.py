# Step 1: Convert base-5 number 3112 to decimal
base_5_number = '3112'
decimal_number = 0
base_5_length = len(base_5_number)

for i in range(base_5_length):
    digit = int(base_5_number[base_5_length - i - 1])
    decimal_number += digit * (5 ** i)

# Step 2: Convert the decimal number to base-7
base_7_number = ''
while decimal_number > 0:
    remainder = decimal_number % 7
    base_7_number = str(remainder) + base_7_number
    decimal_number //= 7

print(base_7_number)