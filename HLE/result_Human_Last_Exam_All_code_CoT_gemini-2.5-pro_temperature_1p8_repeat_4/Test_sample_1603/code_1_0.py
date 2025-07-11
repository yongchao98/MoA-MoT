import string

def solve_cipher():
    """
    Solves a two-step cipher (Caesar then Substitution) to find a quote
    from a Christopher Nolan movie and the character who said it.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    
    # Hint 1: The substitution key
    # The keyword is CHRISTOPENLA (from CHRISTOPHER NOLAN), followed by remaining letters.
    substitution_key = "CHRISTOPENLABDFGJKMQUVWXYZ"
    plain_alphabet = string.ascii_uppercase

    # Hint 2: The number 5 is important, implying a Caesar shift.
    # We will try a backward shift of 5, as is common in puzzles.
    shift = -5

    # --- Step 1: Apply the Caesar Cipher (-5) ---
    caesar_shifted_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            # Calculate the shifted position, handling wraparound
            original_pos = ord(char) - ord('A')
            shifted_pos = (original_pos + shift + 26) % 26
            shifted_char = chr(ord('A') + shifted_pos)
            caesar_shifted_message += shifted_char
        else:
            # Keep spaces as they are
            caesar_shifted_message += char

    # --- Step 2: Apply the Substitution Cipher ---
    # Create the decoding map from the key alphabet to the plain alphabet
    decode_map = str.maketrans(substitution_key, plain_alphabet)
    
    # Translate the result of the Caesar shift using the substitution key
    decoded_message = caesar_shifted_message.translate(decode_map)
    
    # Print the decoded quote
    print(f"The decoded quote is: {decoded_message}")

solve_cipher()
<<<Alfred Pennyworth>>>