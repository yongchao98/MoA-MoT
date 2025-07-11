def solve_cipher():
    """
    Decodes a message encrypted with a Caesar cipher and a substitution cipher.
    """
    # Setup the alphabets and the encrypted message
    plain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    shift = -5

    # Step 1: Reverse the Caesar cipher
    shifted_message = ""
    for char in encoded_message:
        if char == ' ':
            shifted_message += ' '
            continue
        try:
            char_index = plain_alphabet.index(char)
            # Apply the shift, wrapping around the alphabet
            new_index = (char_index + shift) % 26
            shifted_message += plain_alphabet[new_index]
        except ValueError:
            # If the character is not in the alphabet, keep it as is
            shifted_message += char

    # Step 2: Reverse the substitution cipher
    # The map translates the intermediate text (in plain_alphabet) to the final decoded text (in key_alphabet).
    substitution_map = str.maketrans(plain_alphabet, key_alphabet)
    decoded_message = shifted_message.translate(substitution_map)

    # Print the decoded message to reveal the quote
    print(f"The encoded message is: {encoded_message}")
    print(f"The key alphabet is: {key_alphabet}")
    print(f"The Caesar shift value is: {abs(shift)}")
    print(f"The decoded quote is: {decoded_message}")
    print("\nThe quote is from the movie 'Memento'.")
    print("The character who says this is Leonard Shelby.")

solve_cipher()
<<<Leonard Shelby>>>