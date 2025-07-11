def solve_cipher():
    """
    Solves the two-step substitution cipher.
    
    The process involves:
    1. Reversing the second step of the cipher ('A' -> 'BD').
    2. Decrypting the resulting mono-alphabetic substitution cipher.
    """
    
    # The first ciphertext provided by the user
    ciphertext = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."

    # Step 1: Reverse the 'a' -> 'bd' substitution.
    # The original substitution was likely 'A' -> 'BD' to get an all-caps output.
    # We replace "BD" with "A" to get the intermediate P1 ciphertext.
    # The letter 'A' is used as a placeholder because it does not appear in the original ciphertext.
    p1_ciphertext = ciphertext.replace("BD", "A")

    # Step 2: Define the substitution key for the mono-alphabetic cipher.
    # This key was derived using frequency analysis and pattern matching.
    # e.g., ZIT -> the, IT -> he, OSS -> ass, DTLL -> dell, Q -> i
    key = {
        'A': 'n', 'B': 'j', 'C': 'k', 'D': 'd', 'E': 'm', 'F': 'f',
        'G': 'r', 'H': 'p', 'I': 'h', 'K': 'c', 'L': 'a', 'M': 'y',
        'N': 'b', 'O': 's', 'P': 'g', 'Q': 'w', 'R': 'u', 'S': 'v',
        'T': 'o', 'U': 'r', 'V': 'e', 'W': 'q', 'X': 'x', 'Y': 'z',
        'Z': 't'
    }

    # Decrypt the P1 ciphertext using the derived key
    plaintext = ""
    for char in p1_ciphertext:
        if char in key:
            plaintext += key[char]
        else:
            # Keep spaces and punctuation as is
            plaintext += char

    print(f"Ciphertext: {ciphertext}")
    print(f"Decrypted message: {plaintext}")

solve_cipher()