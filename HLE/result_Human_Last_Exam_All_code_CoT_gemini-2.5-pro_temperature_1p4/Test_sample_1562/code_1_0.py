import textwrap

def solve_crypto_phrase():
    """
    Solves the encrypted phrase puzzle by mapping hex codes to letters
    and reconstructing the original sentence.
    """

    # The final mapping discovered through deduction
    code_map = {
        '45a0': 'F',
        'afc5': 'U',
        '3a8d': 'T',
        '0fa7': 'R',
        '529f': 'E',
        '7ba1': 'H',
        'd4c2': 'O',
        '644f': 'L',
        'fb1c': 'D',
        '05d9': 'S',
        'bdb9': 'L',
        '591e': 'S',
        '8f61': 'A',
        '09be': 'M',
        '237f': 'N'
    }

    encrypted_phrase = "45a0afc53a8dafc50fa7529f 7ba1d4c2644ffb1c05d9 bdb9d4c2591e05d9529f05d9 d4c245a0 644fd4c2bdb9237f fb1c529f8f61fb1c fb1c0fa7529f8f6109be05d9"
    encrypted_words = encrypted_phrase.split(' ')
    
    decrypted_words = []

    print("Breaking down the decipherment process word by word:")
    print("-" * 50)

    for i, word in enumerate(encrypted_words):
        # Split the word into 4-character codes
        codes = textwrap.wrap(word, 4)
        
        decrypted_word = ""
        equation_parts = []
        
        for code in codes:
            letter = code_map.get(code, '?')
            decrypted_word += letter
            equation_parts.append(f"{code}({letter})")
        
        decrypted_words.append(decrypted_word)
        print(f"Word {i+1}: {' '.join(equation_parts)}")
        print(f"Result: {decrypted_word}\n")

    final_phrase = ' '.join(decrypted_words)
    print("-" * 50)
    print("The final deciphered phrase is:")
    print(final_phrase)

solve_crypto_phrase()
<<<FUTURE HOLDS LOSSES OF LONE DEAD DREAMS>>>