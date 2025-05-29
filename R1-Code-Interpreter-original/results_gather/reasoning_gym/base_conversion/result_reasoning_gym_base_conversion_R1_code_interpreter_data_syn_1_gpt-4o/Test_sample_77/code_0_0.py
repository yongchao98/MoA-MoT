# Convert base-5 number 2121 to base-10
base_5_number = '2121'
base_5 = 5
base_10_number = sum(int(digit) * (base_5 ** idx) for idx, digit in enumerate(reversed(base_5_number)))

# Convert base-10 number to base-9
base_9 = 9
base_9_number = ''
while base_10_number > 0:
    remainder = base_10_number % base_9
    base_9_number = str(remainder) + base_9_number
    base_10_number //= base_9

print(base_9_number)