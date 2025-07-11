import string

def decode_message():
    """
    Decodes a message that has been encrypted with a substitution cipher
    and a Caesar cipher.
    """
    # The provided encoded message, key, and shift value
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    cipher_key = "CHRISTOPENLABDFGJKMQUVWXYZ"
    caesar_shift = 5

    # Standard alphabet for reference
    alphabet = string.ascii_uppercase

    # Step 1: Reverse the Caesar cipher by shifting back by 5
    # The original encoding process likely involved shifting letters forward by 5,
    # so we shift them back to reverse it.
    caesar_decoded_message = ""
    for char in encoded_message:
        if char in alphabet:
            original_index = alphabet.find(char)
            # Apply the reverse shift, wrapping around the alphabet if necessary
            decoded_index = (original_index - caesar_shift + 26) % 26
            caesar_decoded_message += alphabet[decoded_index]
        else:
            # Keep spaces and other non-alphabetic characters as they are
            caesar_decoded_message += char

    # Step 2: Reverse the substitution cipher
    # The key 'CHRISTOPEN...' maps to the alphabet 'ABCDEF...'. To decode,
    # we create a reverse mapping from the key back to the alphabet.
    # For example, 'C' maps to 'A', 'H' to 'B', and so on.
    decryption_map = str.maketrans(cipher_key, alphabet)
    final_decoded_message = caesar_decoded_message.translate(decryption_map)

    # Output the steps and the final decoded message
    # Let's show the equation part by part as requested.
    print("Decoding process:")
    print(f"Original encoded message: {encoded_message}")

    # Print the Caesar decoding part of the "equation"
    caesar_equation_parts = []
    for char in encoded_message:
        if char in alphabet:
            caesar_equation_parts.append(f"CaesarDecode({char})")
        else:
            caesar_equation_parts.append(f"'{char}'")
    print(f"Step 1 (Caesar Decode): {' + '.join(caesar_equation_parts)} = {caesar_decoded_message}")
    
    # Print the substitution decoding part of the "equation"
    sub_equation_parts = []
    for char in caesar_decoded_message:
         if char in alphabet:
            sub_equation_parts.append(f"SubDecode({char})")
         else:
            sub_equation_parts.append(f"'{char}'")
    print(f"Step 2 (Substitution Decode): {' + '.join(sub_equation_parts)} = {final_decoded_message}")
    
    print("\nFinal decoded quote:")
    print(final_decoded_message)

decode_message()