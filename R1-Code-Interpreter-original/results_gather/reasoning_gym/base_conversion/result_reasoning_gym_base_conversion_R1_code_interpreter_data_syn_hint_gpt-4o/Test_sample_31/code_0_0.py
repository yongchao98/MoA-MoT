# Convert base-9 number 555 to base-10
base_9_number = '555'
base_9 = 9
decimal_number = sum(int(digit) * (base_9 ** idx) for idx, digit in enumerate(reversed(base_9_number)))

# Convert the decimal number to base-6
base_6 = 6
base_6_number = ''
while decimal_number > 0:
    remainder = decimal_number % base_6
    base_6_number = str(remainder) + base_6_number
    decimal_number //= base_6

print(base_6_number)