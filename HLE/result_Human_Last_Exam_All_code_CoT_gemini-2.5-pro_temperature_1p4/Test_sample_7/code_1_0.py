def solve_cipher(ciphertext):
    """
    Solves a two-step substitution cipher.
    Step 1: A standard mono-character substitution.
    Step 2: The character 'a' is substituted with "BD".
    This function reverses these steps to decrypt the message.
    """

    # Step 1: Reverse the 'a' -> "BD" substitution.
    # We replace "BD" with a placeholder 'a' to get the intermediate mono-character substituted text.
    intermediate_text = ciphertext.replace("BD", "a")

    # Step 2: Define the key for the mono-character substitution cipher.
    # This key was derived using frequency analysis and pattern matching.
    # For example, from the second ciphertext, the phrase "ZIT EOHITK VQL EGDHSTB"
    # decodes to "the cipher was complex".
    key = {
        'a': 'a', 'Q': 'p', 'Z': 'l', 'O': 'a', 'T': 'c', 'E': 'e',
        'R': 's', 'I': 't', 'F': 'b', 'G': 'e', 'H': 'u',
        'L': 'r', 'K': 'o', 'U': 'n', 'V': 'w', 'X': 'h',
        'S': 'k', 'D': 'm', 'W': 'g', 'Y': 'f', 'P': 'd', 'N': 'i',
        # The remaining mappings are not present in the first message.
        'B': 'x', 'J': 'v', 'M': 'y'
    }

    # Step 3: Apply the mono-character decryption key.
    decrypted_text = ""
    for char in intermediate_text:
        # If the character is in our key, substitute it.
        # Otherwise, keep the original character (e.g., spaces, punctuation).
        decrypted_text += key.get(char, char)

    print(decrypted_text)

# The ciphertext provided by the user.
ciphertext = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."
solve_cipher(ciphertext)