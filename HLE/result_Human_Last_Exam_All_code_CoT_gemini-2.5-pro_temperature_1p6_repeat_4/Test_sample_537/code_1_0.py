def solve_puzzle():
    """
    This function decodes the word 'anhaeketn' based on a discovered pattern.
    The pattern is that the decoded value of a word is its consonant count
    minus its vowel count.
    """
    word = "anhaeketn"
    vowels = "aeiou"
    
    vowel_count = 0
    consonant_count = 0
    
    for char in word:
        if char in vowels:
            vowel_count += 1
        else:
            consonant_count += 1
            
    decoded_value = consonant_count - vowel_count
    
    # Per the instructions, output each number in the final equation.
    print(f"{consonant_count} - {vowel_count} = {decoded_value}")

solve_puzzle()