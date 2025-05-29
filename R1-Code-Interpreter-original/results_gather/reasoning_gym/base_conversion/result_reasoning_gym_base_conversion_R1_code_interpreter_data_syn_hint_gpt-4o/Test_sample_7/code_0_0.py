# Step 1: Convert base-6 number 554 to decimal
base_6_number = '554'
decimal_number = sum(int(digit) * (6 ** idx) for idx, digit in enumerate(reversed(base_6_number)))

# Step 2: Convert the decimal number to base-8
base_8_number = ''
while decimal_number > 0:
    remainder = decimal_number % 8
    base_8_number = str(remainder) + base_8_number
    decimal_number //= 8

print(base_8_number)