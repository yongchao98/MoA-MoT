# Step 1: Convert base-5 to decimal
base_5_number = '4042'
decimal_number = sum(int(digit) * (5 ** idx) for idx, digit in enumerate(reversed(base_5_number)))

# Step 2: Convert decimal to base-11
base_11_number = ''
while decimal_number > 0:
    remainder = decimal_number % 11
    # Convert remainder to appropriate base-11 digit
    if remainder >= 10:
        base_11_digit = chr(remainder - 10 + ord('a'))
    else:
        base_11_digit = str(remainder)
    base_11_number = base_11_digit + base_11_number
    decimal_number //= 11

print(base_11_number)