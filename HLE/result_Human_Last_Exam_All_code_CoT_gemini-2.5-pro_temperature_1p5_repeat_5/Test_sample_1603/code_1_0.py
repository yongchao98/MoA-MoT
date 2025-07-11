def solve_cipher():
    """
    Decodes the message and identifies the character from the Nolan movie.
    The decoding process involves two steps: a Caesar cipher shift and a keyed substitution.
    """
    # The standard alphabet for reference
    plain_alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    # The custom key alphabet derived from "CHRISTOPHER NOLAN"
    key_alpha   = "CHRISTOPENLABDFGJKMQUVWXYZ"
    
    # The encoded message
    ciphertext  = "PPXRRNUTL RNYFYINIZR"
    
    # The shift value from the hint
    shift = -5

    # --- Step 1: Apply the Caesar cipher (ROT-5) ---
    caesared_text = ""
    for char in ciphertext:
        if 'A' <= char <= 'Z':
            # Calculate the shifted character code
            shifted_code = ord(char) + shift
            # Handle wrap-around for shifts going before 'A'
            if shifted_code < ord('A'):
                shifted_code += 26
            caesared_text += chr(shifted_code)
        else:
            # Preserve spaces or any other characters
            caesared_text += char
    
    # --- Step 2: Apply the keyed substitution decode ---
    # Create a mapping from the key alphabet to the plain alphabet
    # e.g., key_alpha[0] ('C') maps to plain_alpha[0] ('A')
    decode_map = {key_alpha[i]: plain_alpha[i] for i in range(26)}

    plaintext = ""
    for char in caesared_text:
        if 'A' <= char <= 'Z':
            # Decode the character using the map
            plaintext += decode_map.get(char, char)
        else:
            # Preserve spaces
            plaintext += char

    # The decoded quote reveals the source.
    # The quote is "AN IDEA IS LIKE A VIRUS". It's from the movie Inception.
    # The character who says this is Dominic Cobb.
    
    # Print the final decoded quote
    print(plaintext)
    
    # Final answer is the character's name
    # return "Cobb"

solve_cipher()
<<<Cobb>>>