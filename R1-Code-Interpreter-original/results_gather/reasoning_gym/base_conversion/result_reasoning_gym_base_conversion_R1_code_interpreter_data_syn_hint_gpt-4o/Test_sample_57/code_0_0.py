# Step 1: Convert "5d" from base-14 to base-10
base_14_number = "5d"
base_14 = 14
base_10_number = int(base_14_number, base_14)

# Step 2: Convert the base-10 number to base-11
base_11 = 11
digits = "0123456789abcdefghijklmnopqrstuvwxyz"
base_11_number = ""

while base_10_number > 0:
    remainder = base_10_number % base_11
    base_11_number = digits[remainder] + base_11_number
    base_10_number //= base_11

print(base_11_number)