def solve_cipher():
    """
    Solves the two-step cipher puzzle.
    1. Reverses a Caesar cipher with a shift of -5.
    2. Reverses a substitution cipher using the provided key.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    substitution_key = "CHRISTOPENLABDFGJKMQUVWXYZ"
    plain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = -5

    # Step 1: Reverse the Caesar cipher
    print("Step 1: Reversing the Caesar Cipher by shifting each character by -5.")
    print("Equation: new_char_index = (original_char_index - 5) % 26\n")
    
    caesar_decoded_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            original_index = ord(char) - ord('A')
            # The equation for the shift
            new_index = (original_index + shift + 26) % 26
            new_char = chr(ord('A') + new_index)
            print(f"{char}({original_index}) - 5 => {new_char}({new_index})")
            caesar_decoded_message += new_char
        else:
            # For spaces
            print(f"'{char}' => '{char}'")
            caesar_decoded_message += char
    
    print(f"\nResult after Caesar shift: {caesar_decoded_message}\n")

    # Step 2: Reverse the substitution cipher
    print("Step 2: Reversing the substitution cipher.")
    print(f"Using substitution key: {substitution_key}\n")
    
    # Create the decoding map from the key
    decode_map = {key_char: plain_char for key_char, plain_char in zip(substitution_key, plain_alphabet)}

    final_plaintext = ""
    for char in caesar_decoded_message:
        if char in decode_map:
            decoded_char = decode_map[char]
            # Outputting the 'equation' for substitution
            print(f"{char} => {decoded_char}")
            final_plaintext += decoded_char
        else:
            print(f"'{char}' => '{char}'")
            final_plaintext += char
            
    print("\n--- The Decoded Quote is ---")
    print(final_plaintext)

solve_cipher()
<<<Eames>>>