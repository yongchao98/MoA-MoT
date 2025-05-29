# Convert base-6 number 554 to base-10
base_6_number = '554'
base_10_number = sum(int(digit) * (6 ** idx) for idx, digit in enumerate(reversed(base_6_number)))

# Convert base-10 number to base-8
base_8_number = ''
while base_10_number > 0:
    remainder = base_10_number % 8
    base_8_number = str(remainder) + base_8_number
    base_10_number //= 8

print(base_8_number)