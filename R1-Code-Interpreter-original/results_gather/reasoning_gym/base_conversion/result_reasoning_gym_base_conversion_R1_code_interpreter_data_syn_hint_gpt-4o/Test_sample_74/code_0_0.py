# Step 1: Convert from Hexadecimal to Decimal
hex_number = '2d7'
decimal_number = int(hex_number, 16)

# Step 2: Convert from Decimal to Base-9
def decimal_to_base9(n):
    if n == 0:
        return '0'
    digits = []
    while n:
        digits.append(int(n % 9))
        n //= 9
    return ''.join(str(x) for x in digits[::-1])

base9_number = decimal_to_base9(decimal_number)
print(base9_number)