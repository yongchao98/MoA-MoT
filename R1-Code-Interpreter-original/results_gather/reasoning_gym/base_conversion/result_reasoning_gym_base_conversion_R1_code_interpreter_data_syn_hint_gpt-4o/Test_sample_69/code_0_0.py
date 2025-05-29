# Step 1: Convert base-15 number '29a' to base-10
base_15_number = '29a'
base_15 = 15

# Convert to base-10
decimal_value = int(base_15_number, base_15)

# Step 2: Convert the decimal value to base-5
base_5 = 5
base_5_number = ''

while decimal_value > 0:
    remainder = decimal_value % base_5
    base_5_number = str(remainder) + base_5_number
    decimal_value //= base_5

print(base_5_number)