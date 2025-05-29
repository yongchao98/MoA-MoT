# Step 1: Convert base-11 number to decimal
base_11_number = "78"
decimal_number = int(base_11_number, 11)

# Step 2: Convert decimal number to base-4
base_4_number = ""
while decimal_number > 0:
    remainder = decimal_number % 4
    base_4_number = str(remainder) + base_4_number
    decimal_number //= 4

print(base_4_number)