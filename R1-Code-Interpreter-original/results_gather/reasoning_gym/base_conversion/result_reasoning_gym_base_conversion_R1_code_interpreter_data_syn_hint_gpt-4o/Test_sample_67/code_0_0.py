# Convert base-8 number 635 to base-10
base_8_number = '635'
base_10_number = sum(int(digit) * (8 ** idx) for idx, digit in enumerate(reversed(base_8_number)))

# Convert base-10 number to base-5
base_5_number = ''
while base_10_number > 0:
    base_5_number = str(base_10_number % 5) + base_5_number
    base_10_number //= 5

print(base_5_number)