# Convert decimal to base-15
def decimal_to_base15(decimal_number):
    if decimal_number == 0:
        return '0'
    
    digits = []
    while decimal_number > 0:
        remainder = decimal_number % 15
        if remainder >= 10:
            # Convert remainder to corresponding letter for base > 10
            digits.append(chr(remainder - 10 + ord('a')))
        else:
            digits.append(str(remainder))
        decimal_number //= 15
    
    # The digits are in reverse order
    return ''.join(reversed(digits))

base_15_number = decimal_to_base15(982)
print(base_15_number)