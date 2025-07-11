def solve_cipher():
    """
    Solves the two-step substitution cipher.
    Step 1: Replaces "BD" with "a".
    Step 2: Applies the deduced monoalphabetic substitution key.
    """
    ciphertext = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."

    # Step 1: Reverse the 'a' -> 'BD' substitution
    intermediate_text = ciphertext.replace("BD", "a")

    # Step 2: Define the monoalphabetic substitution key found through analysis
    # (Ciphertext -> Plaintext)
    substitution_key = {
        'Q': 'l', 'Z': 'i', 'O': 'f', 'T': 'e', 'E': 'k', 'R': 'n',
        'I': 'h', 'F': 'a', 'G': 'r', 'L': 'g', 'K': 'd', 'U': 'o',
        'V': 'y', 'X': 'u', 'S': 's', 'D': 'c', 'Y': 'p', 'H': 'v',
        'W': 'b', 'N': 'm', 'P': 'w', 'M': 'z', 'C': 'x', 'J': 'j', 'B': 'q'
    }

    plaintext = ""
    for char in intermediate_text:
        # Check if the character is in our substitution key (i.e., an uppercase letter)
        if char in substitution_key:
            plaintext += substitution_key[char]
        else:
            # Keep non-alphabetic characters (spaces, punctuation) and our placeholder 'a' as is
            plaintext += char
            
    print("Deciphered message:")
    print(plaintext)

solve_cipher()