def solve_cipher():
    """
    Decodes a message that has been encrypted with a substitution cipher
    and a Caesar cipher.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    substitution_key = "CHRISTOPENLABDFGJKMQUVWXYZ"
    standard_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift_amount = 5

    # Step 1: Apply the inverse Caesar shift (-5) to the encoded message.
    shifted_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            # Find the character's position (0-25) and shift it back.
            original_pos = ord(char) - ord('A')
            shifted_pos = (original_pos - shift_amount) % 26
            shifted_message += standard_alphabet[shifted_pos]
        else:
            # Keep non-alphabetic characters like spaces as they are.
            shifted_message += char

    # Step 2: Apply the inverse substitution using the key.
    # Create a decryption map from the substitution key to the standard alphabet.
    # e.g., 'C' (from key) maps to 'A', 'H' maps to 'B', etc.
    decryption_map = {substitution_key[i]: standard_alphabet[i] for i in range(26)}

    decoded_message = ""
    for char in shifted_message:
        if 'A' <= char <= 'Z':
            # Look up the character in our map to get the original letter.
            decoded_message += decryption_map[char]
        else:
            decoded_message += char
    
    # Print the decoded quote
    print(f"The decoded quote is: \"{decoded_message}\"")

solve_cipher()
<<<Cobb>>>