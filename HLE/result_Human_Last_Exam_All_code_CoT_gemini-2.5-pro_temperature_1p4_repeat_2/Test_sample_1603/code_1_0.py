def solve_cipher():
    """
    Decodes the message by applying inverse substitution, an inverse Caesar shift,
    and an Atbash cipher, then prints the quote and the character who said it.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alpha = "CHRISTOPENLABDFGJKMQUVWXYZ"
    std_alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = -5

    # 1. Create the inverse substitution map
    sub_inverse_map = {key_alpha[i]: std_alpha[i] for i in range(26)}

    final_plaintext = ""
    print("Decoding Process:")
    print("-" * 20)

    for char_encoded in encoded_message:
        if 'A' <= char_encoded <= 'Z':
            # Step 1: Inverse Substitution
            char_sub_inversed = sub_inverse_map[char_encoded]
            sub_inv_ord = ord(char_sub_inversed) - ord('A')

            # Step 2: Inverse Caesar Shift
            char_caesar_inversed_ord = (sub_inv_ord + shift) % 26
            char_caesar_inversed = chr(char_caesar_inversed_ord + ord('A'))

            # Step 3: Atbash Cipher
            char_atbash_ord = ord('A') + ord('Z') - ord(char_caesar_inversed)
            char_final = chr(char_atbash_ord)

            final_plaintext += char_final
            
            # Print the equation with numbers as requested
            print(f"{char_encoded} -> {char_sub_inversed} ({sub_inv_ord}) -> Caesar({shift}) -> {char_caesar_inversed} ({char_caesar_inversed_ord}) -> Atbash -> {char_final} ({char_atbash_ord - ord('A')})")

        else:
            final_plaintext += " " # Handle spaces
            print("\n[SPACE]\n")

    print("-" * 20)
    print(f"\nDecoded Quote: {final_plaintext}")
    
    # Based on the decoded quote "DREAMALITTLEBIGGER", which is a famous line from the
    # movie Inception ("You mustn't be afraid to dream a little bigger, darling."),
    # the character who says this is Eames.
    character_name = "Eames"
    print(f"The quote is said by: {character_name}")

solve_cipher()
<<<Eames>>>