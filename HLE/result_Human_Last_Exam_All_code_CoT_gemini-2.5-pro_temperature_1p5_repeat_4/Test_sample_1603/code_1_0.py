def solve_cipher():
    """
    Solves the cipher puzzle to find the quote and the character who said it.
    The cipher is a Vigenere cipher variant using a custom alphabet and the keyword 'TENET'.
    """
    std_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    key_alphabet = 'CHRISTOPENLABDFGJKMQUVWXYZ'
    ciphertext = 'PPXRRNUTL RNYFYINIZR'
    keyword = 'TENET'  # The number 5 hints at this 5-letter keyword from a Nolan movie.

    # Remove spaces for easier processing
    ciphertext_no_spaces = ciphertext.replace(" ", "")
    full_key = (keyword * (len(ciphertext_no_spaces) // len(keyword) + 1))[:len(ciphertext_no_spaces)]

    plaintext = ""
    
    print("Decoding Process:")
    print("Ciphertext: {}".format(ciphertext))
    print("Keyword:    {}".format(keyword))
    print("-" * 25)
    print("Cipher | Key | C_idx | K_idx | P_idx | Plain")
    print("------------------------------------------")

    word_lengths = [len(word) for word in ciphertext.split(' ')]
    char_count = 0

    for i in range(len(ciphertext_no_spaces)):
        cipher_char = ciphertext_no_spaces[i]
        key_char = full_key[i]

        # Find the index of the cipher character in the custom key_alphabet
        c_idx = key_alphabet.find(cipher_char)
        # Find the index of the key character in the standard alphabet
        k_idx = std_alphabet.find(key_char)

        if c_idx != -1:
            # Calculate the plaintext index
            p_idx = (c_idx - k_idx + 26) % 26
            plain_char = std_alphabet[p_idx]
            
            # Add a space based on original word lengths
            char_count += 1
            if char_count in [sum(word_lengths[:j+1]) for j in range(len(word_lengths)-1)]:
                 plaintext += plain_char + " "
            else:
                 plaintext += plain_char
                 
            # Print the step-by-step decoding for the final equation
            print("{:^6} | {:^3} | {:^5} | {:^5} | {:^5} | {:^5}".format(
                cipher_char, key_char, c_idx, k_idx, p_idx, plain_char))
        else:
            plaintext += cipher_char
    
    print("\nDecoded Message: " + plaintext)

solve_cipher()
<<<The Protagonist>>>