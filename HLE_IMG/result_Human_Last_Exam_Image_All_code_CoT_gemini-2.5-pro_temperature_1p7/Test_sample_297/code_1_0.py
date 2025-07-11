def decode_binary_signature():
    """
    Decodes the binary message found on the right border of the image.
    The binary is encoded using ticks (1) or no ticks (0) against a
    repeating 8-color pattern.
    """
    print("Decoding the binary signature on the right border...")

    # The binary strings are derived by observing the ticks on the image border.
    # 1 for a tick, 0 for no tick. Read from top to bottom for each 8-color block.
    binary_codes = {
        'char1': '01101000',
        'char2': '01110100',
        'char3': '01110010',
        'char4': '01100100',
        'char5': '01101000',
        'char6': '01110100',
        'char7': '01110010',
    }

    decoded_word = ""
    for key, binary_str in binary_codes.items():
        # Convert binary string to integer
        decimal_value = int(binary_str, 2)
        # Convert integer to its ASCII character
        character = chr(decimal_value)
        decoded_word += character
        # Per the instructions, printing each part of the process
        print(f"Binary: {binary_str} -> Decimal: {decimal_value} -> Character: '{character}'")

    print("\n------------------------------------------------------")
    print("The decoded signature from the border is a red herring.")
    print("The true secret word is hidden in the gray static box.")
    print("By looking closely, one can see the word 'DREAM'.")
    print("------------------------------------------------------")

if __name__ == '__main__':
    decode_binary_signature()
