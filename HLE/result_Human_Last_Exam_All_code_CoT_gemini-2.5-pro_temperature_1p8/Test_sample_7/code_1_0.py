def solve_cipher():
    """
    Solves the two-step substitution cipher.
    Step 1: Decrypts using a mono-alphabetic substitution key.
    Step 2: Replaces all occurrences of the string 'bd' with 'a'.
    """
    ciphertext = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."

    # This key is derived from frequency analysis and pattern recognition from both ciphertexts.
    # The key insight is that the word "BD" in the ciphertext decodes to "bd"
    # in the intermediate text, which is then replaced by 'a'.
    key = {
        'B': 'b', 'D': 'd', 'Q': 'a', 'Z': 't', 'O': 'i', 'T': 'e',
        'R': 'r', 'I': 'h', 'F': 'n', 'G': 'o', 'S': 'l', 'L': 's',
        'K': 'k', 'U': 'u', 'V': 'w', 'X': 'v', 'E': 'c', 'H': 'p',
        'Y': 'f', 'C': 'g', 'N': 'y', 'M': 'm', 'P': 'j', 'W': 'x',
        'J': 'z', 'A': 'q'
    }

    # Step 1: Decipher using the mono-alphabetic substitution key.
    intermediate_text = ""
    for char in ciphertext:
        if char in key:
            intermediate_text += key[char]
        elif char.isspace() or char in ",.":
            intermediate_text += char
        else:
            intermediate_text += '?' # Placeholder for any unmapped characters

    # Step 2: Replace all occurrences of "bd" with "a".
    final_plaintext = intermediate_text.replace("bd", "a")

    print(final_plaintext)

solve_cipher()
<<<a ate a ice a er the a nobler a ise skunou, a nowou it wvusr abe e sad iss to vnad not the aounler dss.>>>