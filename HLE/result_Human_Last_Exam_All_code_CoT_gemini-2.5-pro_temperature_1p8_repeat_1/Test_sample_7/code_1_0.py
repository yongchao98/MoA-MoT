def solve_and_print_decryption():
    """
    Solves a two-step substitution cipher.
    The first step is a mono-character substitution.
    The second step replaces 'a' with 'BD'.
    This function reverses the process to find the original plaintext.
    """
    
    # The first ciphertext that needs to be deciphered
    cipher_text = "BD QZOT BD OEBD TR ZIT BD FGZZTR BD OZT LZKOFU, BD FGVOFU OZ VGXSR ZQBD T LBD OSS ZG XFBD FGZ ZIT ZQFUSTR DTLL."
    
    # Step 1: Reverse the substitution of 'a' with 'BD'
    intermediate_text = cipher_text.replace("BD", "a")

    # Step 2: Define the monoalphabetic substitution key
    # This key was derived by analyzing both ciphertexts. Strong clues included:
    # - ZIT -> the
    # - EOHITK -> cipher
    # - ZG -> to
    # - OZ -> it
    # - DTZIGR -> method
    # - WTTF HQLLTR RGVF -> been passed down
    # - VGXSR -> would
    # - ZQFUSTR DTLL -> tangled mess
    key = {
        'A': 'v', 'B': 'x', 'C': 'j', 'D': 'm', 'E': 'c', 'F': 'n', 'G': 'o', 
        'H': 'p', 'I': 'h', 'J': 'k', 'K': 'r', 'L': 's', 'M': 'z', 'N': 'y',
        'O': 'i', 'P': 'q', 'Q': 'a', 'R': 'd', 'S': 'l', 'T': 'e', 'U': 'g', 
        'V': 'w', 'W': 'b', 'X': 'u', 'Y': 'f', 'Z': 't'
    }

    deciphered_text = ""
    # Loop through the intermediate text and substitute characters based on the key
    for char in intermediate_text:
        if 'A' <= char <= 'Z':
            deciphered_text += key.get(char, "?")
        else:
            deciphered_text += char
            
    print(deciphered_text)

solve_and_print_decryption()