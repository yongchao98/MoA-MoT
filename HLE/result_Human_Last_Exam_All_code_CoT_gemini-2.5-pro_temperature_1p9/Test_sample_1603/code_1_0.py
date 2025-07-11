def decode_nolan_quote():
    """
    Decodes a message encrypted with a Keyed Caesar Cipher.
    The process is:
    1. Reverse the Caesar shift on the ciphertext.
    2. Reverse the substitution using the custom key alphabet.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    std_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = 5

    # Create the reverse substitution map (from key alphabet to standard alphabet)
    decryption_map = {key_alphabet[i]: std_alphabet[i] for i in range(len(key_alphabet))}

    print(f"Decoding message: \"{encoded_message}\"")
    print(f"Using key alphabet: \"{key_alphabet}\" and a reverse shift of {shift}.")
    print("-" * 30)

    decoded_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            # Step 1: Reverse Caesar Shift
            # Find the position of the character in the standard alphabet (0-25)
            char_pos = ord(char) - ord('A')
            # Shift backwards, wrapping around the alphabet if necessary
            shifted_pos = (char_pos - shift + 26) % 26
            shifted_char = std_alphabet[shifted_pos]

            # Step 2: Reverse Substitution
            decoded_char = decryption_map[shifted_char]
            
            # Print the step-by-step logic for each character
            print(f"'{char}' -> shift back {shift} -> '{shifted_char}' -> substitute -> '{decoded_char}'")
            decoded_message += decoded_char
        else:
            # Handle spaces or other characters
            print("\nHandling space\n")
            decoded_message += char

    print("-" * 30)
    print(f"The fully decoded message is: \"{decoded_message}\"")

    # The quote is "YOURE WAITING FOR A TRAIN", said by Mal in Inception.
    character_name = "Mal"
    print(f"The quote is from the movie Inception, spoken by the character: {character_name}")

decode_nolan_quote()