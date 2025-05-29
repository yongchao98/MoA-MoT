# Convert base-14 number '41c' to decimal
base_14_number = '41c'
decimal_value = int(base_14_number, 14)

# Convert the decimal value to binary
binary_value = bin(decimal_value)[2:]

print(binary_value)