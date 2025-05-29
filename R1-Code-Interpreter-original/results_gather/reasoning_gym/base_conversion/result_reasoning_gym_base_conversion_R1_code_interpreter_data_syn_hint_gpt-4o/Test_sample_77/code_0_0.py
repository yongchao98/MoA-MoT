# Convert base-5 number 2121 to base-10
base_5_number = "2121"
base_10_number = sum(int(digit) * (5 ** idx) for idx, digit in enumerate(reversed(base_5_number)))

# Convert base-10 number to base-9
base_9_number = ""
while base_10_number > 0:
    remainder = base_10_number % 9
    base_9_number = str(remainder) + base_9_number
    base_10_number //= 9

print(base_9_number)