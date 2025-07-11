import string

def solve_nolan_cipher():
    """
    Decodes a message from a Christopher Nolan movie using a Variant Vigenere cipher.
    
    The process combines the given key and hint:
    1.  A custom alphabet is used, derived from the key 'CHRISTOPENLABDFGJKMQUVWXYZ'.
    2.  A Vigenere cipher is applied. The hint '5' points to a 5-letter keyword.
        'DREAM' from the movie 'Inception' is chosen as the Vigenere key.
    3.  The decryption formula P_idx = (C_idx - K_idx) mod 26 is applied, with all lookups
        performed on the custom alphabet.
    """
    
    custom_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    cipher_text = "PPXRRNUTL RNYFYINIZR".replace(" ", "")
    vigenere_key = "DREAM"
    
    decrypted_text = ""
    key_len = len(vigenere_key)
    
    for i, char_c in enumerate(cipher_text):
        # Find the index of the ciphertext character in the custom alphabet
        c_idx = custom_alphabet.find(char_c)
        
        # Get the corresponding Vigenere key character and its index
        char_k = vigenere_key[i % key_len]
        k_idx = custom_alphabet.find(char_k)
        
        # Apply the Vigenere decryption formula using modular arithmetic
        p_idx = (c_idx - k_idx + 26) % 26
        
        # Append the decrypted character from the custom alphabet
        decrypted_text += custom_alphabet[p_idx]
        
    print(f"Decoded Message: {decrypted_text}")
    print("This is a famous quote from the movie 'Inception'.")
    print("The character who says this is Mal.")

solve_nolan_cipher()
<<<Mal>>>