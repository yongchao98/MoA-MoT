def binary_to_hex(binary_str):
    # Convert binary to integer
    decimal_value = int(binary_str, 2)
    # Convert integer to hexadecimal
    hex_value = hex(decimal_value)[2:]  # [2:] to remove the '0x' prefix
    return hex_value

binary_number = "100010010"
hexadecimal_number = binary_to_hex(binary_number)
print(hexadecimal_number)