def base13_to_base10(number):
    base13_digits = '0123456789abc'
    base10 = 0
    for i, digit in enumerate(reversed(number)):
        base10 += base13_digits.index(digit) * (13 ** i)
    return base10

def base10_to_base5(number):
    if number == 0:
        return '0'
    base5 = ''
    while number > 0:
        base5 = str(number % 5) + base5
        number //= 5
    return base5

# Convert base-13 'b6' to base-10
base10_number = base13_to_base10('b6')

# Convert base-10 number to base-5
base5_number = base10_to_base5(base10_number)

print(base5_number)