# Convert base-6 number 3442 to base-10
base_6_number = "3442"
base_10_number = sum(int(digit) * (6 ** idx) for idx, digit in enumerate(reversed(base_6_number)))

# Convert base-10 number to base-7
base_7_number = ""
while base_10_number > 0:
    remainder = base_10_number % 7
    base_7_number = str(remainder) + base_7_number
    base_10_number //= 7

print(base_7_number)