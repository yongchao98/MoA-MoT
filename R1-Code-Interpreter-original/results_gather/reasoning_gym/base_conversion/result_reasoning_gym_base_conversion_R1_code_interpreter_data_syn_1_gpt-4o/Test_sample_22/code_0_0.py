# Step 1: Convert base-4 to base-10
base_4_number = '32213'
base_10_number = 0
base_4 = 4

# Calculate base-10 equivalent
for index, digit in enumerate(reversed(base_4_number)):
    base_10_number += int(digit) * (base_4 ** index)

# Step 2: Convert base-10 to base-9
base_9_number = ''
base_9 = 9

# Calculate base-9 equivalent
while base_10_number > 0:
    remainder = base_10_number % base_9
    base_9_number = str(remainder) + base_9_number
    base_10_number //= base_9

print(base_9_number)