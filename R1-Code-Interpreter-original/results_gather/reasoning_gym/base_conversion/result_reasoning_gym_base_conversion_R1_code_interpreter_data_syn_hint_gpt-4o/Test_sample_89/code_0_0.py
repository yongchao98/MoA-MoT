# Convert base-5 number 324 to decimal
base_5_number = '324'
decimal_number = sum(int(digit) * (5 ** idx) for idx, digit in enumerate(reversed(base_5_number)))

# Convert decimal number to base-11
base_11_number = ''
while decimal_number > 0:
    remainder = decimal_number % 11
    if remainder >= 10:
        # Convert remainder to corresponding letter for bases greater than 10
        base_11_number = chr(remainder - 10 + ord('a')) + base_11_number
    else:
        base_11_number = str(remainder) + base_11_number
    decimal_number //= 11

print(base_11_number)