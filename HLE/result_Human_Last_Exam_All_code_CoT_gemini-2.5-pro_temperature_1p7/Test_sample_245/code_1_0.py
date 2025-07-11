import re
from russtress import Stresser

def find_words_with_non_initial_stress():
    """
    Analyzes a Russian text to find multi-syllable words where the stress
    is not on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    
    # Initialize the stresser. Note: the first time this runs, it may download a model.
    try:
        stresser = Stresser()
    except Exception as e:
        print(f"Could not initialize the Stresser. Please ensure 'russtress' is installed ('pip install russtress').")
        print(f"Error: {e}")
        return

    # A string containing all Russian vowels.
    vowels = "аоуыэяёюие"
    
    # Extract words, preserving dashes within words like 'кто-то'.
    original_words = re.findall(r'\b[\w-]+\b', text)
    
    # This list will hold the final words that meet the criteria.
    result_words = []
    
    # We iterate through the original words to maintain order and original casing.
    for original_word in original_words:
        # Processing is done on the lowercase version of the word.
        word_lower = original_word.lower()
        
        # 1. Count syllables. Words with 1 or 0 syllables are excluded.
        syllable_count = sum(1 for char in word_lower if char in vowels)
        if syllable_count <= 1:
            continue
            
        # 2. Find the index of the first vowel.
        first_vowel_index = -1
        for i, char in enumerate(word_lower):
            if char in vowels:
                first_vowel_index = i
                break
        
        # 3. Determine the index of the stressed vowel.
        stressed_vowel_index = -1
        # The letter 'ё' is always stressed, which simplifies finding the stress.
        if 'ё' in word_lower:
            stressed_vowel_index = word_lower.find('ё')
        else:
            # For other words, use the russtress library.
            stressed_word = stresser.stress_word(word_lower)
            # The library adds a combining acute accent (U+0301) after the stressed vowel.
            stress_accent = "\u0301"
            if stress_accent in stressed_word:
                # The index of the stressed vowel is right before the accent mark.
                stressed_vowel_index = stressed_word.find(stress_accent) - 1

        # 4. Compare the index of the stressed vowel to the first vowel.
        # If they are different, the stress is not on the first syllable.
        # We also check that we successfully found a stress.
        if stressed_vowel_index != -1 and stressed_vowel_index != first_vowel_index:
            result_words.append(original_word)

    # Print the final list, comma-separated.
    print(", ".join(result_words))

find_words_with_non_initial_stress()
<<<Шашлык, горелым, прибежал, сосед, дошёл, его>>>