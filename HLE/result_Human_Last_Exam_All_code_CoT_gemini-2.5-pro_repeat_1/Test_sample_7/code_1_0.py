def solve_cipher():
    """
    Solves the two-step substitution cipher.
    Step 1: Reverses the 'a' -> 'BD' substitution.
    Step 2: Applies the deduced monoalphabetic substitution key to decrypt the message.
    """
    ciphertext = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."

    # Step 1: Reverse the 'a' -> 'BD' substitution. Let's use 'a' as the placeholder.
    intermediate_text = ciphertext.replace("BD", "a")

    # Step 2: Define the monoalphabetic substitution key.
    # This key was derived through frequency analysis and pattern matching.
    # 'a' represents the letter substituted by 'BD'.
    key = {
        'a': 'i', 'B': 'x', 'C': 'p', 'D': 'm', 'E': 'c', 'F': 'n',
        'G': 'o', 'H': 'p', 'I': 'h', 'K': 'r', 'L': 's', 'N': 'y',
        'O': 'u', 'P': 'd', 'Q': 'a', 'R': 'd', 'S': 'l', 'T': 'e',
        'U': 'g', 'V': 'w', 'W': 'b', 'X': 'u', 'Y': 'f', 'Z': 't'
    }

    # Decrypt the intermediate text to get the final plaintext.
    plaintext = ""
    for char in intermediate_text:
        if 'A' <= char <= 'Z':
            plaintext += key.get(char, '?')
        elif char == 'a':
            plaintext += key.get('a', '?')
        else:
            plaintext += char

    print(plaintext)

solve_cipher()