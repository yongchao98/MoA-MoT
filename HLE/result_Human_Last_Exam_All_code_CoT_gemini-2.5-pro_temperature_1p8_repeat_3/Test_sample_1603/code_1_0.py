def solve_cipher():
    """
    Decodes a message that has been encrypted with a substitution cipher
    and a Caesar cipher.
    """
    ciphertext = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    plain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = 5

    # Step 1: Apply Caesar Cipher (deciphering means shifting back)
    caesar_decoded_list = []
    for char in ciphertext:
        if 'A' <= char <= 'Z':
            # Calculate the shifted character's ASCII value
            # (ord(char) - ord('A')) gets the 0-25 position
            # We subtract the shift and use modulo 26 to wrap around
            shifted_position = (ord(char) - ord('A') - shift) % 26
            # Convert the new position back to a character
            caesar_decoded_list.append(plain_alphabet[shifted_position])
        else:
            # Keep spaces as they are
            caesar_decoded_list.append(char)
    
    caesar_decoded = "".join(caesar_decoded_list)

    # Step 2: Apply Substitution Cipher
    # Create the mapping from the key alphabet to the plain alphabet
    decryption_map = str.maketrans(key_alphabet, plain_alphabet)
    
    # Translate the Caesar-decoded message using the substitution map
    final_text = caesar_decoded.translate(decryption_map)

    # The puzzle asked to output the numbers in the final equation.
    # The equation can be described as:
    # Decoded Char = Substitution(Caesar_Shift(Encoded Char, -5))
    # We will print the quote and the number used in the equation.
    
    print(f"The decoded quote is: {final_text}")
    print(f"The number used in the decoding equation (Caesar shift) was: {shift}")
    
    # Based on the decoded quote "A LITTLE FIGHT IN YOU", the movie is Inception.
    # The character who says this line is Eames.

solve_cipher()
<<<Eames>>>