def convert_base_10_to_base_n(number, base):
    if number == 0:
        return "0"
    
    digits = []
    while number:
        remainder = number % base
        if remainder >= 10:
            # Convert remainder to corresponding letter for bases > 10
            digits.append(chr(remainder - 10 + ord('a')))
        else:
            digits.append(str(remainder))
        number //= base
    
    # The digits are collected in reverse order
    return ''.join(reversed(digits))

# Convert 331 from base-10 to base-11
result = convert_base_10_to_base_n(331, 11)
print(result)