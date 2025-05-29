def convert_to_base(n, base):
    digits = "0123456789abcdefghijklmnopqrstuvwxyz"
    if n == 0:
        return "0"
    result = []
    while n > 0:
        remainder = n % base
        result.append(digits[remainder])
        n = n // base
    return ''.join(reversed(result))

# Convert 778 from base-10 to base-13
base_10_number = 778
base_13_number = convert_to_base(base_10_number, 13)
print(base_13_number)