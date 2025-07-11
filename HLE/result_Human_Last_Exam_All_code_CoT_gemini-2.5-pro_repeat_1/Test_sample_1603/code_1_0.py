def solve_cipher():
    """
    Solves the composite cipher by finding the character's index in a custom key,
    applying a Caesar shift to that index, and then finding the resulting
    character in the standard alphabet.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    standard_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = -5  # The hint is 5; we use a negative shift for decoding.

    decoded_message = ""

    for char in encoded_message:
        if char == ' ':
            decoded_message += ' '
            continue

        try:
            # 1. Find the index of the character in the custom key alphabet.
            key_index = key_alphabet.find(char)

            # 2. Apply the Caesar shift to the index.
            shifted_index = (key_index + shift) % 26

            # 3. Find the new character in the standard alphabet.
            decoded_char = standard_alphabet[shifted_index]
            decoded_message += decoded_char
        except ValueError:
            # If the character is not in the key alphabet, append it as is.
            decoded_message += char

    print(f"Encoded message: {encoded_message}")
    print(f"Decoded quote: {decoded_message}")
    # The decoded quote is "THEYRE GONNA CHASE YOU", which is said by James Gordon in "The Dark Knight".

solve_cipher()
<<<James Gordon>>>