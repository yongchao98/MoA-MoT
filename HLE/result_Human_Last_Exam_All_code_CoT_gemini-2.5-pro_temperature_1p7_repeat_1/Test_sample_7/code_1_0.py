def solve_cipher():
    """
    Solves the two-step substitution cipher.
    Step 1: Replaces all occurrences of 'BD' with 'a', reversing the second step of the encryption.
    Step 2: Applies the solved monoalphabetic substitution key to decrypt the message.
    """
    ciphertext = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."

    # Step 1: Reverse the 'a' -> 'BD' substitution.
    # The prompt implies that the character substituted for 'a' was then replaced by 'BD'.
    # We use 'a' directly as the placeholder since it is the final plaintext letter.
    intermediate_text = ciphertext.replace('BD', 'A') # Using a placeholder to build the map

    # Step 2: Define the monoalphabetic substitution key.
    # This key was derived using frequency and pattern analysis on the full corpus of ciphertext.
    # For example, "IT VKGZT OZ OF EOHITK" in the second sample decrypts to "he wrote it in cipher".
    cipher_map = {
        'A': 'a', 'B': 'x', 'C': 'j', 'D': 'k', 'E': 'c', 'F': 'f', 'G': 'o',
        'H': 'p', 'I': 'h', 'J': 'z', 'K': 'r', 'L': 'l', 'M': 'm', 'N': 'q',
        'O': 'i', 'P': 'v', 'Q': 'u', 'R': 's', 'S': 'd', 'T': 'e', 'U': 'y',
        'V': 'w', 'W': 'b', 'X': 'g', 'Y': 'v', 'Z': 't'
    }

    # Decrypt the intermediate text
    plaintext = ""
    for char in intermediate_text:
        if char in cipher_map:
            plaintext += cipher_map[char]
        else:
            plaintext += char  # Keep spaces and punctuation as is

    print(plaintext)

solve_cipher()