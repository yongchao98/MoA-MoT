def base14_to_base10(number):
    base14_digits = '0123456789abcdef'
    base10 = 0
    for i, digit in enumerate(reversed(number)):
        base10 += base14_digits.index(digit) * (14 ** i)
    return base10

def base10_to_base9(number):
    if number == 0:
        return '0'
    base9_digits = []
    while number > 0:
        base9_digits.append(str(number % 9))
        number //= 9
    return ''.join(reversed(base9_digits))

# Convert base-14 number 'a8' to base-10
base10_number = base14_to_base10('a8')

# Convert the base-10 number to base-9
base9_number = base10_to_base9(base10_number)

print(base9_number)