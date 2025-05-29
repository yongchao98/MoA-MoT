def decimal_to_base5(n):
    if n == 0:
        return '0'
    digits = []
    while n:
        digits.append(int(n % 5))
        n //= 5
    return ''.join(str(x) for x in digits[::-1])

# Convert hexadecimal 4f to decimal
hex_number = '4f'
decimal_number = int(hex_number, 16)

# Convert decimal to base-5
base5_number = decimal_to_base5(decimal_number)
print(base5_number)