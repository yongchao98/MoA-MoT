def convert_base_10_to_base_n(number, base):
    if number == 0:
        return "0"
    
    digits = []
    while number > 0:
        remainder = number % base
        if remainder >= 10:
            # Convert remainder to corresponding letter for bases > 10
            digits.append(chr(remainder - 10 + ord('a')))
        else:
            digits.append(str(remainder))
        number = number // base
    
    # Reverse the digits to get the correct base-n representation
    return ''.join(reversed(digits))

# Convert 778 from base-10 to base-13
result = convert_base_10_to_base_n(778, 13)
print(result)