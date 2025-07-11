# First, you need to install the 'russtress' library.
# You can do this by running the following command in your terminal or command prompt:
# pip install russtress

import re
from russtress import Stress

def solve():
    """
    Finds and prints words from a given Russian text that have more than one syllable
    and are not stressed on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    
    # Russian vowels for syllable counting.
    vowels = "аеёиоуыэюя"
    
    # Initialize the stress detector.
    stress = Stress()
    
    # Use a regex to split the text by spaces and punctuation, keeping words intact.
    # This handles commas, spaces, and hyphens as delimiters.
    words = re.split(r'[\s,\-]+', text)
    
    result_words = []
    
    for word in words:
        # Clean any leftover punctuation, e.g., a period at the end of the sentence.
        clean_word = word.strip('.,-')
        
        # Skip empty strings that might result from splitting.
        if not clean_word:
            continue
            
        lower_word = clean_word.lower()
        
        # 1. Filter out one-syllable words by counting vowels.
        vowel_count = sum(1 for char in lower_word if char in vowels)
        if vowel_count <= 1:
            continue
            
        # 2. Determine the stress position.
        first_vowel_index = -1
        for i, char in enumerate(lower_word):
            if char in vowels:
                first_vowel_index = i
                break
        
        # This should not happen for multi-syllable words, but as a safeguard.
        if first_vowel_index == -1:
            continue
            
        stressed_vowel_index = -1
        # The letter 'ё' is always stressed.
        if 'ё' in lower_word:
            stressed_vowel_index = lower_word.find('ё')
        else:
            # Use the library to find the stress.
            accentuated_word = stress.accentuate(lower_word)
            # The accent mark '́' is a combining character, placed after the vowel.
            accent_pos = accentuated_word.find('́')
            if accent_pos != -1:
                # The index of the vowel in the original word is the position before the accent mark.
                stressed_vowel_index = accent_pos - 1

        # 3. Check if the stress is NOT on the first syllable.
        # This means the index of the stressed vowel is not the same as the index of the first vowel.
        if stressed_vowel_index != -1 and stressed_vowel_index != first_vowel_index:
            result_words.append(clean_word)
            
    # Print the final comma-separated list.
    print(", ".join(result_words))

solve()
<<<Шашлык, горелым, прибежал, сосед, дошёл, его>>>