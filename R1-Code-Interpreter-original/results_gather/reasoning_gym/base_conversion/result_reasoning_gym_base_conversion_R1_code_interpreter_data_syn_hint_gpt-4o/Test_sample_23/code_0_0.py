# Convert base-5 number 332 to base-10
base_5_number = '332'
base_5 = 5
base_10_number = sum(int(digit) * (base_5 ** idx) for idx, digit in enumerate(reversed(base_5_number)))

# Convert base-10 number to base-4
base_4 = 4
base_4_number = ''
while base_10_number > 0:
    remainder = base_10_number % base_4
    base_4_number = str(remainder) + base_4_number
    base_10_number //= base_4

print(base_4_number)