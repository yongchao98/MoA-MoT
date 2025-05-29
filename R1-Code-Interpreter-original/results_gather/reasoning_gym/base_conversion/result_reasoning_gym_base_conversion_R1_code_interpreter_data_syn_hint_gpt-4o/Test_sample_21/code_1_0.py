# Step 1: Convert base-15 number '2c1' to decimal
base_15_number = '2c1'
decimal_value = int(base_15_number, 15)

# Step 2: Convert the decimal value to hexadecimal
hexadecimal_value = hex(decimal_value)[2:]  # [2:] to remove the '0x' prefix

print(hexadecimal_value)