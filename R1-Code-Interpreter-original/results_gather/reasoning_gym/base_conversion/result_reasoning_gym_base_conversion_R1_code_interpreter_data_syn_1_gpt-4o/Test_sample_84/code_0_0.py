def base13_to_decimal(base13_num):
    base13_digits = '0123456789abc'
    decimal_value = 0
    for i, digit in enumerate(reversed(base13_num)):
        decimal_value += base13_digits.index(digit) * (13 ** i)
    return decimal_value

def decimal_to_base11(decimal_num):
    base11_digits = '0123456789a'
    if decimal_num == 0:
        return '0'
    base11_num = ''
    while decimal_num > 0:
        remainder = decimal_num % 11
        base11_num = base11_digits[remainder] + base11_num
        decimal_num //= 11
    return base11_num

# Convert base-13 number "4c9" to decimal
decimal_value = base13_to_decimal('4c9')

# Convert the decimal value to base-11
base11_value = decimal_to_base11(decimal_value)

print(base11_value)