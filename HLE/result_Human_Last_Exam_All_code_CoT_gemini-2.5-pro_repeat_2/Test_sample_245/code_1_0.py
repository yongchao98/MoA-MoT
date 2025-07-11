# The following script requires the 'udar' library.
# You can install it by running: pip install udar

import re
import udar

def find_first_vowel_index(word, vowels):
    """Finds the index of the first vowel in a word."""
    for i, char in enumerate(word.lower()):
        if char in vowels:
            return i
    return -1

def solve():
    """
    Finds and prints Russian words from a given text that have more than one syllable
    and are not stressed on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    vowels = "аеёиоуыэюя"
    
    # Use regex to find all words containing Russian letters, preserving their original form.
    words = re.findall(r'\b[а-яА-ЯёЁ-]+\b', text)
    
    result_words = []
    
    for word in words:
        word_lower = word.lower()
        
        # Count vowels to determine if the word is polysyllabic.
        num_vowels = sum(1 for char in word_lower if char in vowels)
        
        # Exclude one-syllable words.
        if num_vowels <= 1:
            continue
            
        stress_index = -1
        
        # The Cyrillic letter 'ё' is always stressed.
        if 'ё' in word_lower:
            stress_index = word_lower.find('ё')
        else:
            # Get the stressed version of the word using the udar library.
            # We take the first variant, which is typically the most common.
            stressed_variants = udar.accent(word_lower)
            if stressed_variants:
                stressed_word = stressed_variants[0]
                # Find the position of the stress mark (a combining acute accent).
                stress_mark_pos = stressed_word.find('́')
                if stress_mark_pos != -1:
                    # The stressed vowel is the character right before the accent mark.
                    stress_index = stress_mark_pos - 1

        # If we couldn't determine the stress, skip the word.
        if stress_index == -1:
            continue
            
        # Find the index of the very first vowel in the word.
        first_vowel_index = find_first_vowel_index(word, vowels)
        
        # If the index of the stressed vowel is greater than the index of the first vowel,
        # it means the stress is not on the first syllable.
        if stress_index > first_vowel_index:
            result_words.append(word)
            
    print(", ".join(result_words))

solve()
<<<Шашлык, горелым, прибежал, сосед, дошёл, его>>>