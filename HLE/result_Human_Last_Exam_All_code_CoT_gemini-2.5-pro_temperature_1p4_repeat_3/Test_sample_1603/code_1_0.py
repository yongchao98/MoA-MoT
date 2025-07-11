import string

def solve_cipher():
    """
    Solves a multi-step cipher to find a quote from a Christopher Nolan movie.
    The decryption process involves reversing a substitution cipher and a Caesar cipher.
    """
    
    # 1. Define the knowns from the problem description.
    ciphertext = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    standard_alphabet = string.ascii_uppercase
    shift_key = 5

    # 2. Create the decryption map for the substitution cipher.
    # This maps each character in the key_alphabet back to its standard alphabet equivalent.
    # For example, the first letter of the key 'C' maps to 'A', 'H' maps to 'B', and so on.
    substitution_decryption_map = {key_alphabet[i]: standard_alphabet[i] for i in range(26)}

    # 3. Decrypt the message by first reversing the substitution, then reversing the Caesar cipher.
    # This corresponds to an encryption process of Caesar shift then substitution.
    
    # First, apply the substitution decryption.
    substituted_text = ""
    for char in ciphertext:
        if char in substitution_decryption_map:
            substituted_text += substitution_decryption_map[char]
        else:
            substituted_text += char  # Preserve spaces

    # Second, apply the Caesar shift decryption.
    plaintext = ""
    for char in substituted_text:
        if char in standard_alphabet:
            original_index = standard_alphabet.find(char)
            # Shift back by 5, wrapping around the alphabet if necessary.
            decrypted_index = (original_index - shift_key) % 26
            plaintext += standard_alphabet[decrypted_index]
        else:
            plaintext += char # Preserve spaces
    
    # 4. Print the final decrypted quote.
    print(f"The decoded quote is: {plaintext}")

solve_cipher()
<<<Cobb>>>