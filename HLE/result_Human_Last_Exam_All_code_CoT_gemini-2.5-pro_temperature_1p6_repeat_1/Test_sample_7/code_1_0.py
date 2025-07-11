def solve_cipher():
    """
    Solves the two-step substitution cipher.
    Step 1: Reverses the 'a' -> 'BD' substitution.
    Step 2: Applies the derived key for the monoalphabetic substitution.
    """

    target_ciphertext = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."

    # The substitution key was derived using frequency analysis and pattern recognition
    # on the full ciphertext provided in the prompt, after reversing the 'BD'->'a' step.
    # Key Derivations:
    # ZIT -> the => Z:t, I:h, T:e
    # VQL -> was => V:w, Q:a, L:s
    # RGVF -> down => R:d, G:o, F:n
    # WTTF -> been => W:b
    # YQZITK -> father => Y:f, K:r
    # LNLZTD -> system => N:y, D:m
    # HQLLTR -> passed => H:p
    # EGDHSTB -> complex => E:c, S:l, B:x
    # VGXSR -> would => X:u
    # OY -> if => O:i
    # DTLLQUT -> message => U:g
    # EGRTWKTQ(A) -> combreaka (from intermediate text), suggests 'A' in the ciphertext maps to 'k'.
    
    key = {
        'A': 'k', 'B': 'x', 'C': 'v', 'D': 'm', 'E': 'c', 'F': 'n', 'G': 'o', 'H': 'p',
        'I': 'h', 'K': 'r', 'L': 's', 'M': 'z', 'N': 'y', 'O': 'i', 'P': 'j', 'Q': 'a',
        'R': 'd', 'S': 'l', 'T': 'e', 'U': 'g', 'V': 'w', 'W': 'b', 'X': 'u', 'Y': 'f', 'Z': 't'
    }

    decoded_parts = []
    # Split the original ciphertext by space to process word by word
    for word in target_ciphertext.split(' '):
        # Isolate punctuation
        punctuation = ''
        if not word.isalpha() and word:
            last_char = word[-1]
            if last_char in [',', '.']:
                punctuation = last_char
                word = word[:-1]

        # First, reverse the special 'BD' -> 'a' substitution
        temp_word = word.replace('BD', 'a')
        
        # Second, apply the monoalphabetic substitution key
        decoded_word = ""
        for char in temp_word:
            # If the character is 'a', it's already decoded.
            if char == 'a':
                decoded_word += 'a'
            # Otherwise, use the key to find the plaintext character.
            elif char.upper() in key:
                decoded_word += key[char.upper()]
            # Keep any unmapped characters to identify them.
            else:
                decoded_word += char

        decoded_parts.append(decoded_word + punctuation)

    final_message = " ".join(decoded_parts)
    print(final_message)

solve_cipher()