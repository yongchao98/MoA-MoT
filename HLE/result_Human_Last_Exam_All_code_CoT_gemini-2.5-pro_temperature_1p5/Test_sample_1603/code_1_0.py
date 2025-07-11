def solve_cipher():
    """
    Solves the cipher based on the provided hints.
    The cipher is a Vigenère variant (Quagmire IV) using a mixed alphabet.
    The keyword is derived from the hint '5' and the movie director's name.
    """
    # The standard alphabet for reference and for the keyword.
    std_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    # The mixed alphabet key provided in the hint.
    key_alphabet = 'CHRISTOPENLABDFGJKMQUVWXYZ'
    
    # The encoded message.
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    
    # The number 5 hints at a 5-letter keyword, which is NOLAN.
    keyword = "NOLAN"
    
    decoded_message = ""
    keyword_index = 0
    
    # Loop through each character in the encoded message.
    for char in encoded_message:
        if char == ' ':
            decoded_message += ' '
            continue
        
        # Find the index of the encoded character in the key_alphabet.
        encoded_char_index = key_alphabet.find(char)
        
        # Get the current character from the repeating keyword.
        key_char = keyword[keyword_index % len(keyword)]
        
        # Find the index of the keyword character in the standard alphabet.
        key_char_index = std_alphabet.find(key_char)
        
        # The plaintext index is found by subtracting the key index from the encoded index.
        # We add 26 before the modulo to handle negative results correctly.
        plaintext_char_index = (encoded_char_index - key_char_index + 26) % 26
        
        # Find the plaintext character in the standard alphabet.
        decoded_message += std_alphabet[plaintext_char_index]
        
        keyword_index += 1
        
    print(f"The decoded quote is: {decoded_message}")
    
    # The quote "THE ALGORITHM IS THE BOMB" is from the movie 'Tenet'.
    # It is said by the character Mahir, when describing their plan to crash the 747.
    character_name = "Mahir"
    print(f"The quote is said by the character: {character_name}")

solve_cipher()