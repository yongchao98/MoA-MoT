# Step 1: Convert base-8 number 426 to base-10
base_8_number = '426'
base_10_number = int(base_8_number, 8)

# Step 2: Convert the base-10 number to base-9
def convert_to_base_n(number, base):
    if number == 0:
        return '0'
    digits = []
    while number:
        digits.append(int(number % base))
        number //= base
    return ''.join(str(x) for x in digits[::-1])

base_9_number = convert_to_base_n(base_10_number, 9)
print(base_9_number)