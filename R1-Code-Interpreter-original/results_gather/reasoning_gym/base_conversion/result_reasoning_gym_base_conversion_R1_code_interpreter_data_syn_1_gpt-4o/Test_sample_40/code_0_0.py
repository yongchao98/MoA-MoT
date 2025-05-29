def convert_base_10_to_base_n(number, base):
    if number == 0:
        return "0"
    digits = []
    while number:
        digits.append(int(number % base))
        number //= base
    return ''.join(str(x) for x in digits[::-1])

# Base-10 number
base_10_number = 75

# Convert to base-11
base_11_number = convert_base_10_to_base_n(base_10_number, 11)
print(base_11_number)