def binary_to_hexadecimal(binary_str):
    # Convert binary to integer
    decimal_value = int(binary_str, 2)
    # Convert integer to hexadecimal
    hexadecimal_value = hex(decimal_value)[2:]  # [2:] to remove the '0x' prefix
    return hexadecimal_value

binary_number = "100010010"
hexadecimal_number = binary_to_hexadecimal(binary_number)
print(hexadecimal_number)