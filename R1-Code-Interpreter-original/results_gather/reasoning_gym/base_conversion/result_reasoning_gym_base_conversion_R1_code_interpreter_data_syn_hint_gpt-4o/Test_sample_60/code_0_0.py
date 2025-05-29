def convert_base_10_to_base_n(number, base):
    if number == 0:
        return "0"
    
    digits = []
    while number:
        digits.append(int(number % base))
        number //= base
    
    # Convert digits to appropriate base-n representation
    base_n_number = ''.join(str(d) if d < 10 else chr(d - 10 + ord('a')) for d in reversed(digits))
    return base_n_number

# Convert 431 from base-10 to base-11
result = convert_base_10_to_base_n(431, 11)
print(result)