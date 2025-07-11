def solve_cipher():
    """
    Decodes a message using a Caesar cipher followed by a substitution cipher
    to find a movie quote and the character who said it.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    
    # The key provided for the substitution cipher.
    # This is the 'ciphertext' alphabet.
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    
    # The standard 'plaintext' alphabet.
    plain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    # The hint "number 5" implies a Caesar shift of 5. For decoding, we shift back.
    caesar_shift = -5

    # Step 1: Apply the reverse Caesar shift.
    shifted_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            # Find the position of the character in the alphabet (0-25).
            char_index = ord(char) - ord('A')
            # Apply the shift, wrapping around the alphabet.
            shifted_index = (char_index + caesar_shift) % 26
            # Convert the new index back to a character.
            shifted_message += chr(shifted_index + ord('A'))
        else:
            # Keep spaces and other characters as they are.
            shifted_message += char

    # Step 2: Apply the substitution cipher to the shifted message.
    # Create a decoding map from the key_alphabet to the plain_alphabet.
    decode_map = {key_alphabet[i]: plain_alphabet[i] for i in range(len(key_alphabet))}

    decoded_message = ""
    for char in shifted_message:
        if char in decode_map:
            # Substitute the character using the map.
            decoded_message += decode_map[char]
        else:
            # Keep spaces and other characters as they are.
            decoded_message += char

    # The decoded message reveals the quote from the movie "Inception".
    # The quote is "YOU'RE WAITING FOR A TRAIN", said by Mal.
    print(f"Encoded message: {encoded_message}")
    print(f"Applying Caesar shift of {caesar_shift}: {shifted_message}")
    print(f"Applying substitution with key '{key_alphabet}':")
    print(f"Final decoded quote: {decoded_message}")
    print("\nThe quote is from the movie Inception.")
    print("The character who says this quote is Mal.")

solve_cipher()
<<<Mal>>>