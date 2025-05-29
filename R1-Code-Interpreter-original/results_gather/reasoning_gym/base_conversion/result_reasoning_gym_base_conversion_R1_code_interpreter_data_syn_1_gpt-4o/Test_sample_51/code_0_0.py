def base15_to_decimal(base15_num):
    base15_digits = '0123456789abcdef'
    decimal_value = 0
    for i, digit in enumerate(reversed(base15_num)):
        decimal_value += base15_digits.index(digit) * (15 ** i)
    return decimal_value

def decimal_to_base9(decimal_num):
    if decimal_num == 0:
        return '0'
    base9_digits = []
    while decimal_num > 0:
        base9_digits.append(str(decimal_num % 9))
        decimal_num //= 9
    return ''.join(reversed(base9_digits))

# Convert base-15 number "23d" to decimal
decimal_value = base15_to_decimal('23d')

# Convert the decimal value to base-9
base9_value = decimal_to_base9(decimal_value)

print(base9_value)