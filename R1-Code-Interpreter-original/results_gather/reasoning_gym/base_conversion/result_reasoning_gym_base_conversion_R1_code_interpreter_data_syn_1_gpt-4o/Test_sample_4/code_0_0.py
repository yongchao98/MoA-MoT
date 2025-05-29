# Step 1: Convert binary to decimal
binary_number = '101010000'
decimal_number = int(binary_number, 2)

# Step 2: Convert decimal to base-14
def decimal_to_base14(n):
    if n == 0:
        return '0'
    digits = []
    while n:
        remainder = n % 14
        if remainder >= 10:
            # Convert remainder to corresponding letter for base > 10
            digits.append(chr(remainder - 10 + ord('a')))
        else:
            digits.append(str(remainder))
        n //= 14
    return ''.join(reversed(digits))

base14_number = decimal_to_base14(decimal_number)
print(base14_number)