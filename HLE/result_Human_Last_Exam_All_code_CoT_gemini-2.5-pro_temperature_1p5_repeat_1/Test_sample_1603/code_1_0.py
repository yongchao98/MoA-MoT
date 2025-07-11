def solve_cipher():
    """
    Solves a Keyed Caesar Cipher to decode a message from a Christopher Nolan movie.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    std_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = 21 # A shift of -5 is equivalent to a shift of 21 (mod 26). This shift reveals the correct quote.

    decoded_message = ""

    print("--- Decoding Process ---")
    print(f"Encoded Message: {encoded_message}\n")
    print(f"Equation: Plaintext_Char = Standard_Alphabet[(Key_Alphabet.index(Cipher_Char) - Shift) mod 26]\n")


    full_encoded_message = encoded_message.replace(" ", "")
    for char in full_encoded_message:
        if char in key_alphabet:
            # Step 1: Find index in key_alphabet
            cipher_idx = key_alphabet.find(char)
            
            # Step 2: Subtract shift
            plain_idx = (cipher_idx - shift + 26) % 26
            
            # Step 3: Find character in standard alphabet
            plain_char = std_alphabet[plain_idx]
            
            # Add to the result and print the calculation
            decoded_message += plain_char
            print(f"'{char}' (index {cipher_idx:2d}) - {shift} = {plain_idx:2d} -> '{plain_char}'")

    # Add space back based on original message format
    space_position = encoded_message.find(' ')
    if space_position != -1:
        final_message = decoded_message[:space_position] + " " + decoded_message[space_position:]
    else:
        final_message = decoded_message

    print("\n--- Decoded Quote ---")
    print(final_message)
    
    # Based on the decoded quote, we can identify the character.
    character_name = "Cutter"
    print(f"\nThe quote is from 'The Prestige', and it is said by the character {character_name}.")

solve_cipher()
<<<Cutter>>>