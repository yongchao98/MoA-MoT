def base15_to_decimal(base15_num):
    base15_digits = '0123456789abcdef'
    decimal_value = 0
    for i, digit in enumerate(reversed(base15_num)):
        decimal_value += base15_digits.index(digit) * (15 ** i)
    return decimal_value

def decimal_to_base5(decimal_num):
    if decimal_num == 0:
        return '0'
    base5_digits = []
    while decimal_num > 0:
        base5_digits.append(str(decimal_num % 5))
        decimal_num //= 5
    return ''.join(reversed(base5_digits))

# Convert base-15 number "29a" to decimal
decimal_value = base15_to_decimal('29a')

# Convert the decimal value to base-5
base5_value = decimal_to_base5(decimal_value)

print(base5_value)