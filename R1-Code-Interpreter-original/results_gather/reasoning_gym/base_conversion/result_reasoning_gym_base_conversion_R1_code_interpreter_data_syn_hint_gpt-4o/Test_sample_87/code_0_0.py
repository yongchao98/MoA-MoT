# Convert base-6 number 4220 to base-10
base_6_number = "4220"
base_10_number = 0
base_6 = 6

# Convert from base-6 to base-10
for index, digit in enumerate(reversed(base_6_number)):
    base_10_number += int(digit) * (base_6 ** index)

# Convert from base-10 to base-4
base_4_number = ""
base_4 = 4

while base_10_number > 0:
    remainder = base_10_number % base_4
    base_4_number = str(remainder) + base_4_number
    base_10_number //= base_4

print(base_4_number)