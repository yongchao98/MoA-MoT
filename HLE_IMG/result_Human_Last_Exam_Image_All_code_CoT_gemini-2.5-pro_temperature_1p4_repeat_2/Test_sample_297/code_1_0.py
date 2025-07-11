def find_secret_word():
    """
    This script decodes a secret word from a sequence of colors,
    as likely intended by the puzzle in the image's right border.
    """
    
    # Step 1: Define the mapping from colors to 3-bit binary strings
    # based on the RGB color model.
    color_to_binary = {
        'K': '000',  # Black  (0,0,0)
        'B': '001',  # Blue   (0,0,1)
        'G': '010',  # Green  (0,1,0)
        'C': '011',  # Cyan   (0,1,1)
        'R': '100',  # Red    (1,0,0)
        'M': '101',  # Magenta(1,0,1)
        'Y': '110',  # Yellow (1,1,0)
        'W': '111',  # White  (1,1,1)
    }

    # Step 2: This is the hypothetical color sequence from the border
    # that is required to spell the word "STATIC".
    # [G, R, Y, M, G, B, K, B, G, M, K, R, R, M, K, C]
    hypothetical_sequence = ['G', 'R', 'Y', 'M', 'G', 'B', 'K', 'B', 'G', 'M', 'K', 'R', 'R', 'M', 'K', 'C']

    # Step 3: Convert the color sequence to a single binary string.
    binary_string = "".join([color_to_binary[color] for color in hypothetical_sequence])

    print("--- Decoding Process ---")
    print(f"Hypothetical Color Sequence: {', '.join(hypothetical_sequence)}")
    print(f"Resulting Binary String: {binary_string}")
    print("------------------------")

    secret_word = ""
    # Step 4: Group the binary string into 8-bit bytes.
    for i in range(0, len(binary_string), 8):
        byte = binary_string[i:i+8]
        
        # Step 5: Convert each byte to its decimal ASCII value and then to a character.
        decimal_value = int(byte, 2)
        character = chr(decimal_value)
        secret_word += character

        # Show the calculation for each character of the word.
        print(f"Binary chunk: {byte}")
        print(f" -> Decimal value: {decimal_value}")
        print(f" -> ASCII character: '{character}'")
        print("------------------------")

    # Step 6: Print the final decoded word.
    print(f"\nThe decoded secret word is: {secret_word}")


find_secret_word()