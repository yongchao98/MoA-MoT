def solve_cipher():
    """
    Solves the two-step cipher to find the character's name.
    """
    # Given information
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    standard_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift_amount = 5

    # --- Step 1: Caesar Cipher Decode ---
    # We apply an inverse shift of 5 to the encoded message.
    caesared_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            # Shift character back by 5, wrapping around the alphabet
            shifted_ord = ord(char) - shift_amount
            if shifted_ord < ord('A'):
                shifted_ord += 26
            caesared_message += chr(shifted_ord)
        else:
            caesared_message += char  # Preserve spaces

    # --- Step 2: Substitution Decode ---
    # We map letters from the key_alphabet back to the standard_alphabet.
    decoding_map = dict(zip(key_alphabet, standard_alphabet))

    final_message = ""
    for char in caesared_message:
        if char in decoding_map:
            final_message += decoding_map[char]
        else:
            final_message += char

    # Print the decoded message to reveal the quote
    print(f"Encoded message: {encoded_message}")
    print(f"After Caesar Decode (shift -{shift_amount}): {caesared_message}")
    print(f"After Substitution Decode: {final_message}")
    print("\nThe quote is from the movie 'The Dark Knight'.")
    print("The character who says this quote is Harvey Dent.")

solve_cipher()
<<<Harvey Dent>>>