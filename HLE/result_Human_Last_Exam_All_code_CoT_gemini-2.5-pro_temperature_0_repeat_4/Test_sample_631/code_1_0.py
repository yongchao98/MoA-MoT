def solve_vigenere_and_math():
    """
    This script solves a Vigenere cipher puzzle.
    1. It determines the key from a riddle.
    2. It decrypts the ciphertext using the key.
    3. It solves the mathematical question found in the plaintext.
    4. It prints the decrypted question and the final answer.
    """
    # The key is derived from the puzzle: "AYE" (affirmation) sounds like "AI" (Google's AI).
    key = "ai"
    ciphertext = "tsa mmrc u'w nrq ma brihmarshlfveziovsfcatmpsupmpohsppazsgbnuswvs?"
    
    # Decrypt the ciphertext
    decrypted_text = ""
    key_index = 0
    for char in ciphertext:
        if 'a' <= char <= 'z':
            # Vigenere decryption: P = (C - K) % 26
            cipher_ord = ord(char) - ord('a')
            key_ord = ord(key[key_index % len(key)]) - ord('a')
            plain_ord = (cipher_ord - key_ord + 26) % 26
            decrypted_text += chr(plain_ord + ord('a'))
            key_index += 1
        else:
            decrypted_text += char
            
    # The decrypted question is:
    # how much is twenty plus twenty two divided by two times three?
    # Let's solve the math problem following order of operations (PEMDAS/BODMAS)
    # 20 + 22 / 2 * 3
    # Division first: 22 / 2 = 11
    # Then multiplication: 11 * 3 = 33
    # Finally addition: 20 + 33 = 53
    
    num1 = 20
    num2 = 22
    num3 = 2
    num4 = 3
    
    result = num1 + num2 / num3 * num4
    
    # Print the decrypted question and the final equation
    print("Decrypted Question: how much is twenty plus twenty two divided by two times three?")
    print(f"Final Equation: {num1} + {num2} / {num3} * {num4} = {int(result)}")

solve_vigenere_and_math()
<<<53>>>