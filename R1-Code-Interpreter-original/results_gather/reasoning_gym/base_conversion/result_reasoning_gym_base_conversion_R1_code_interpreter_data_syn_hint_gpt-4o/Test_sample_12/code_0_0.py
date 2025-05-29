def convert_to_base(n, base):
    digits = []
    while n > 0:
        remainder = n % base
        if remainder >= 10:
            # Convert remainder to corresponding letter for bases > 10
            remainder = chr(remainder - 10 + ord('a'))
        digits.append(str(remainder))
        n = n // base
    return ''.join(digits[::-1])

# Convert 418 to base-13
result = convert_to_base(418, 13)
print(result)