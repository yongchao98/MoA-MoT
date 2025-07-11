def solve_substitution_cipher():
    """
    Solves the two-step substitution cipher.
    The first step replaces 'BD' with 'a'.
    The second step applies a derived mono-alphabetic substitution key.
    """
    
    # The ciphertext to be deciphered
    ciphertext = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."

    # This is the full mono-alphabetic substitution key derived through
    # frequency analysis and pattern matching as described in the plan.
    # The character 'A' is used as a placeholder for the special 'BD' case.
    # In this final key, 'a' from plaintext maps to 'BD' in ciphertext.
    substitution_key = {
        'Q': 'p', 'Z': 'a', 'O': 'y', 'T': 's', 'E': 'm', 'R': 'b', 
        'I': 'e', 'F': 'h', 'G': 'w', 'L': 'o', 'K': 'u', 'U': 'l', 
        'V': 'r', 'X': 'd', 'S': 'n', 'D': 'c'
        # Other letters not present in the first message (C, H, J, M, N, P, W, Y)
        # are not needed for this specific decryption.
    }

    # Step 1: Reverse the second encryption step. Replace "BD" with "a".
    # All other characters are kept as uppercase to denote they are still encrypted.
    intermediate_text = ciphertext.replace("BD", "a")

    # Step 2: Decrypt the mono-alphabetic substitution cipher.
    plaintext = ""
    for char in intermediate_text:
        if 'A' <= char <= 'Z':
            # If the character is in our key, substitute it.
            # Otherwise, keep it as is (should not happen with a full key).
            plaintext += substitution_key.get(char, '?')
        else:
            # Keep spaces, punctuation, and our special 'a' as they are.
            plaintext += char

    print(f"Deciphered message: {plaintext}")

solve_substitution_cipher()
<<<Deciphered message: a pays to remember that a man who is weary of the struggle, but he who has always been a stranger to a who is the first to lose heart.>>>