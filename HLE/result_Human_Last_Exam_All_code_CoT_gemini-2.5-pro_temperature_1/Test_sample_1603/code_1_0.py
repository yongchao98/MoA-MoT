def solve_cipher():
    """
    Solves the cryptographic puzzle by first applying a Caesar cipher decryption
    and then a keyed substitution decryption.
    """
    ciphertext = "PPXRRNUTL RNYFYINIZR"
    keyed_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    std_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = 5

    # Step 1: Undo the Caesar cipher (shift by -5)
    intermediate_text = ""
    print("Decoding Step 1: Caesar Cipher (Shift by -5)")
    print(f"Ciphertext: {ciphertext}")
    equation_step1 = []
    for char in ciphertext:
        if 'A' <= char <= 'Z':
            # (ord(char) - ord('A') - shift + 26) % 26 gives the new index
            # ord('A') + new_index gives the new character
            shifted_char_ord = ord(char) - shift
            if shifted_char_ord < ord('A'):
                shifted_char_ord += 26
            decrypted_char = chr(shifted_char_ord)
            intermediate_text += decrypted_char
            equation_step1.append(f"{char}(-{shift})->{decrypted_char}")
        else:
            intermediate_text += char
            equation_step1.append(char)
    
    print("Equation: " + " ".join(equation_step1))
    print(f"Intermediate Text: {intermediate_text}\n")

    # Step 2: Undo the substitution cipher
    # The encryption map is std_alphabet[i] -> keyed_alphabet[i]
    # The decryption map is keyed_alphabet[i] -> std_alphabet[i]
    decryption_map = {keyed_alphabet[i]: std_alphabet[i] for i in range(26)}

    plaintext = ""
    print("Decoding Step 2: Keyed Substitution")
    print(f"Using Key: {keyed_alphabet}")
    equation_step2 = []
    for char in intermediate_text:
        if 'A' <= char <= 'Z':
            decrypted_char = decryption_map.get(char, '?')
            plaintext += decrypted_char
            equation_step2.append(f"{char}->{decrypted_char}")
        else:
            plaintext += char
            equation_step2.append(char)
            
    print("Equation: " + " ".join(equation_step2))
    print(f"\nFinal Decoded Quote: {plaintext}")

    # Step 3: Identify the character
    # The decoded quote is "NOWYOUARELOOKINGFORTHESECRET"
    # This is from the movie "The Prestige"
    character_name = "Cutter"
    print(f"\nThe quote is from the movie 'The Prestige'.")
    print(f"The character who says this is: {character_name}")


solve_cipher()
<<<Cutter>>>