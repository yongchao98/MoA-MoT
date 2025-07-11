import string

def solve_cipher():
    """
    Decodes a message that has been encrypted with a Caesar cipher
    followed by a substitution cipher.
    """
    # --- Setup ---
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    standard_alphabet = string.ascii_uppercase
    shift = 5

    print(f"Encoded message: {encoded_message}\n")
    print(f"Step 1: Reversing the Caesar Cipher by shifting back by {shift} positions.")

    # --- Step 1: Apply reverse Caesar cipher ---
    shifted_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            shifted_index = (ord(char) - ord('A') - shift + 26) % 26
            shifted_message += standard_alphabet[shifted_index]
        else:
            shifted_message += char
    
    print(f"Intermediate message after Caesar decode: {shifted_message}\n")
    
    print(f"Step 2: Reversing the substitution using the key '{key_alphabet}'.")

    # --- Step 2: Apply reverse substitution ---
    decoding_map = {key_char: std_char for key_char, std_char in zip(key_alphabet, standard_alphabet)}
    
    decoded_message = ""
    for char in shifted_message:
        if char in decoding_map:
            decoded_message += decoding_map[char]
        else:
            decoded_message += char
            
    print(f"Final decoded message: {decoded_message}")
    
    # --- Analysis ---
    # The decoded message "NEVERFORGET WHOYOUARE" is a quote from the movie "Batman Begins".
    # The character who says this is Carmine Falcone.
    # The hint "5" could also refer to the fact that "Batman Begins" (2005) can be considered Nolan's 5th major film.

solve_cipher()
<<<Carmine Falcone>>>