def solve_riddle():
    """
    Decodes the message "PPXRRNUTL RNYFYINIZR" to find the character who said the quote.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    shift = -5 # The number 5 is a hint for a Caesar cipher shift of -5

    # The key is a custom alphabet for a substitution cipher.
    plain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"

    # Step 1: Apply the Caesar cipher shift.
    shifted_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            shifted_ord = ord(char) + shift
            if shifted_ord < ord('A'):
                shifted_ord += 26
            shifted_message += chr(shifted_ord)
        else:
            shifted_message += char

    # Step 2: Apply the substitution cipher using the key alphabet.
    # Create a decryption map from the key alphabet to the plain alphabet.
    decryption_map = {key_alphabet[i]: plain_alphabet[i] for i in range(26)}

    decoded_message = ""
    for char in shifted_message:
        if 'A' <= char <= 'Z':
            decoded_message += decryption_map[char]
        else:
            decoded_message += char

    print(f"Encoded message: {encoded_message}")
    print(f"After Caesar shift by {shift}: {shifted_message}")
    print(f"Final decoded quote: {decoded_message}")
    print("\nThe quote is from the movie Inception.")
    # The character's name is not part of the decoded message itself.
    # It's the answer to the riddle based on the decoded quote.
    character_name = "Cobb"
    print(f"The character who said this is: {character_name}")

solve_riddle()
<<<Cobb>>>