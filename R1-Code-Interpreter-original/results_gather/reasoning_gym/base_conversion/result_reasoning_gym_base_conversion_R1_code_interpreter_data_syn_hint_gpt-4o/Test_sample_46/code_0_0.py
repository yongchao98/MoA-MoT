# Step 1: Convert base-5 number 11314 to decimal
base_5_number = '11314'
decimal_number = 0
base_5_length = len(base_5_number)

for i in range(base_5_length):
    digit = int(base_5_number[i])
    power = base_5_length - i - 1
    decimal_number += digit * (5 ** power)

# Step 2: Convert decimal number to base-3
base_3_number = ''
while decimal_number > 0:
    remainder = decimal_number % 3
    base_3_number = str(remainder) + base_3_number
    decimal_number //= 3

print(base_3_number)