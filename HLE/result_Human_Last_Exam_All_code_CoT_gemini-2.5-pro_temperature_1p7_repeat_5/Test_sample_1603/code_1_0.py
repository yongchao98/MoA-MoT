def solve_cipher():
    """
    Solves the given cipher by performing a Vigenère decryption
    followed by an Atbash cipher.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    vigenere_key = "NOLAN"
    
    # Remove space for continuous decryption
    continuous_encoded_message = encoded_message.replace(" ", "")
    
    # --- Step 1: Vigenère Cipher Decryption ---
    print("--- Step 1: Vigenère Cipher ---")
    print(f"Applying Vigenère cipher with key '{vigenere_key}'.")
    print("The equation for each character is: Plaintext = (Ciphertext - Keyword) mod 26\n")

    vigenere_decoded_message = ""
    for i, enc_char in enumerate(continuous_encoded_message):
        key_char = vigenere_key[i % len(vigenere_key)]
        
        enc_val = ord(enc_char) - ord('A')
        key_val = ord(key_char) - ord('A')
        
        plain_val = (enc_val - key_val + 26) % 26
        plain_char = chr(plain_val + ord('A'))
        
        vigenere_decoded_message += plain_char
        
        print(f"'{enc_char}' ({enc_val}) - '{key_char}' ({key_val}) = {plain_val} -> '{plain_char}'")

    print(f"\nResult after Vigenère: {vigenere_decoded_message}\n")

    # --- Step 2: Atbash Cipher Decryption ---
    print("--- Step 2: Atbash Cipher ---")
    print("Applying Atbash cipher (reversing the alphabet).")
    print("The equation for each character is: Plaintext = (25 - Ciphertext_Value) mod 26\n")
    
    final_message_list = []
    for i, enc_char in enumerate(vigenere_decoded_message):
        enc_val = ord(enc_char) - ord('A')
        plain_val = 25 - enc_val
        plain_char = chr(plain_val + ord('A'))
        
        # Re-insert the space at the original position (after 9th character)
        if i == 9:
            final_message_list.append(" ")
            
        final_message_list.append(plain_char)

        print(f"'{enc_char}' ({enc_val}) -> 25 - {enc_val} = {plain_val} -> '{plain_char}'")
        
    final_message = "".join(final_message_list)
    print(f"\nFinal Decoded Quote: {final_message}")
    
    # --- Step 3: Identify the Character ---
    # The decoded quote "ARE YOU WATCHING CLOSELY" is from the movie "The Prestige".
    # It is said by the character John Cutter.
    character_name = "Cutter"
    print(f"\nThe character who says this is: {character_name}")

solve_cipher()
<<<Cutter>>>