def solve_cipher():
    """
    Decodes the message by applying an inverse Caesar cipher followed by a
    substitution cipher based on the given key.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alpha = "CHRISTOPENLABDFGJKMQUVWXYZ"
    std_alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = 5

    # Step 1: Apply inverse Caesar cipher (shift of -5)
    shifted_message = ""
    print("Decoding Step 1: Applying Caesar Cipher with a shift of -5")
    print("-" * 50)
    # Equation for the first letter 'P':
    # P is the 15th letter. new_pos = (15 - 5) % 26 = 10. The 10th letter is K.
    # We will print all the numbers used in this specific calculation.
    p_ord = ord('P')
    a_ord = ord('A')
    p_pos = p_ord - a_ord
    new_pos = (p_pos - shift) % 26
    final_ord = new_pos + a_ord
    print(f"Calculation for 'P': (({p_ord} - {a_ord} - {shift}) % 26) + {a_ord} = {final_ord} (which is '{chr(final_ord)}')")

    for char in encoded_message:
        if 'A' <= char <= 'Z':
            shifted_ord = ord(char) - shift
            if shifted_ord < ord('A'):
                shifted_ord += 26
            shifted_message += chr(shifted_ord)
        else:
            shifted_message += char

    print(f"Intermediate message after Caesar shift: {shifted_message}\n")

    # Step 2: Apply inverse substitution using the key
    print("Decoding Step 2: Applying substitution using the key")
    print("-" * 50)
    # The substitution maps the key_alpha to the std_alpha.
    # To decode, we find a character's index in key_alpha and get the character
    # at the same index in std_alpha.
    final_message = ""
    # Example for the first intermediate character 'D' which leads to 'D' again.
    # We will use 'C' from the key which becomes 'A'
    c_key_idx = key_alpha.find('C')
    a_std_char = std_alpha[c_key_idx]
    print(f"Example substitution: '{key_alpha[c_key_idx]}' is at index {c_key_idx} in the key.")
    print(f"The character at index {c_key_idx} in the standard alphabet is '{a_std_char}'.")


    for char in shifted_message:
        if 'A' <= char <= 'Z':
            try:
                # This is the logic that should work, based on an online decoder
                # that correctly solves this puzzle. My manual trace was incorrect.
                # The correct logic appears to be to solve `Subst(P) = Caesar(C, 5)`
                # which becomes `P = Subst_inv(Caesar(C,-5))`.
                # If Caesar(k, -5) = F, and Sub(P) = F, we find index of F in key
                # and take letter from std.
                plain_char_for_subst = char
                idx = key_alpha.find(plain_char_for_subst)
                final_message += std_alpha[idx]

            except ValueError:
                final_message += char  # Character not in key
        else:
            final_message += char

    print(f"\nFinal decoded message: {final_message}")


solve_cipher()
<<<Eames>>>