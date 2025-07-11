def solve_vigenere_and_riddle():
    """
    Solves the Vigenere cipher and the subsequent riddle.
    """
    ciphertext = "tsa mmrc u'w nrq ma brihmarshlfveziovsfcatmpsupmpohsppazsgbnuswvs?"
    key = "MEENA"
    
    decrypted_text = ""
    key_index = 0
    
    # Decrypt the ciphertext
    for char in ciphertext:
        if 'a' <= char <= 'z':
            shift = ord(key[key_index % len(key)].lower()) - ord('a')
            decrypted_char_code = (ord(char) - ord('a') - shift + 26) % 26
            decrypted_text += chr(decrypted_char_code + ord('a'))
            key_index += 1
        elif 'A' <= char <= 'Z':
            shift = ord(key[key_index % len(key)].upper()) - ord('A')
            decrypted_char_code = (ord(char) - ord('A') - shift + 26) % 26
            decrypted_text += chr(decrypted_char_code + ord('A'))
            key_index += 1
        else:
            decrypted_text += char
            
    print(f"Decrypted question: {decrypted_text}")
    
    # Answer the decrypted question
    president_name = "Richard Milhous Nixon"
    # Count both lower 'i' and upper 'I'
    count_of_i = president_name.lower().count('i')
    
    print(f"The 37th President's full name is: {president_name}")
    print(f"The number of 'i's in the name is: {count_of_i}")

solve_vigenere_and_riddle()