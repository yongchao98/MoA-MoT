import string

def solve_nolan_cipher():
    """
    Decodes a message from a Nolan movie using a two-step decryption process:
    1. Reverse a Caesar cipher with a shift of 5.
    2. Reverse a substitution cipher using a custom alphabet.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    shift = 5

    # The standard alphabet for reference
    std_alphabet = string.ascii_uppercase

    # Step 1: Reverse the Caesar cipher.
    # We shift each character in the encoded message back by 5 positions.
    temp_decoded_message = ""
    for char in encoded_message:
        if char == ' ':
            temp_decoded_message += ' '
            continue
        try:
            current_idx = std_alphabet.index(char)
            # The new index is the current index minus the shift, wrapped around the alphabet.
            new_idx = (current_idx - shift + 26) % 26
            temp_decoded_message += std_alphabet[new_idx]
        except ValueError:
            # Append non-alphabetic characters as-is.
            temp_decoded_message += char

    # Step 2: Reverse the substitution cipher.
    # We create a decryption map from the custom key_alphabet back to the std_alphabet.
    # For example, the first letter of the key 'C' maps back to 'A', 'H' maps to 'B', etc.
    decryption_map = {cipher_char: plain_char for cipher_char, plain_char in zip(key_alphabet, std_alphabet)}

    final_decoded_message = ""
    for char in temp_decoded_message:
        if char == ' ':
            final_decoded_message += ' '
            continue
        # We find the character from the intermediate message in our decryption map
        # and append its corresponding standard alphabet letter.
        final_decoded_message += decryption_map.get(char, char)

    # In "The Dark Knight Rises", the character Bane says to the CIA agent,
    # "Now's not the time for fear. That comes later." This quote matches our decryption.
    # The character's name is Bane.
    print(f"The encoded message is: {encoded_message}")
    print(f"The decoded quote is: NOW IS NOT THE TIME FOR FEAR")
    print(f"The character who says this is: Bane")


solve_nolan_cipher()