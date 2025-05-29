# Step 1: Convert base-11 number '78' to base-10
base_11_number = '78'
base_11 = 11

# Convert to base-10
base_10_number = int(base_11_number, base_11)

# Step 2: Convert base-10 number to base-4
base_4 = 4
base_4_number = ''

while base_10_number > 0:
    remainder = base_10_number % base_4
    base_4_number = str(remainder) + base_4_number
    base_10_number //= base_4

print(base_4_number)