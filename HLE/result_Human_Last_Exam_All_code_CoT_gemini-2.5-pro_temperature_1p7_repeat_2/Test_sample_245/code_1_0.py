# The following script requires the 'russian-stresses' library.
# You can install it by running: pip install russian-stresses
import re
try:
    from russian_stresses import Stresses
except ImportError:
    print("Error: The 'russian-stresses' library is not installed.")
    print("Please install it using: pip install russian-stresses")
    exit()

def find_words_with_non_first_syllable_stress():
    """
    Analyzes a Russian sentence to find words that are not stressed on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    vowels = "аеёиоуыэюя"
    
    # Initialize the stress analysis tool
    stresses = Stresses()
    
    # Use regex to find all Russian words in the text, preserving their order
    words = re.findall(r'\b[а-яА-ЯёЁ-]+\b', text)
    
    result_words = []
    
    for word in words:
        lower_word = word.lower()
        
        # 1. Filter out one-syllable words
        vowel_count = sum(1 for char in lower_word if char in vowels)
        if vowel_count <= 1:
            continue
        
        # 2. Find the index of the first vowel in the word
        first_vowel_index = -1
        for i, char in enumerate(lower_word):
            if char in vowels:
                first_vowel_index = i
                break
        
        # This check is for safety, should not be triggered for words with vowels
        if first_vowel_index == -1:
            continue

        # 3. Determine if the word's stress falls on the first syllable
        is_stressed_on_first = False
        
        # The letter 'ё' is always stressed, handle it as a special case
        if 'ё' in lower_word:
            stressed_vowel_index = lower_word.find('ё')
            if stressed_vowel_index == first_vowel_index:
                is_stressed_on_first = True
        else:
            # For other words, get stress patterns from the library
            stressed_options = stresses.stressify(lower_word)
            
            # Policy: if any possible stress pattern is on the first syllable, exclude the word.
            for option in stressed_options:
                # The stressed vowel is the one just before the accent mark ('\u0301')
                stressed_vowel_idx_in_original = option.find('\u0301') - 1
                if stressed_vowel_idx_in_original == first_vowel_index:
                    is_stressed_on_first = True
                    break # Found a first-syllable stress, no need to check other patterns

        # 4. If after all checks, stress is not on the first syllable, add original word to results
        if not is_stressed_on_first:
            result_words.append(word)
            
    # 5. Print the final list
    print(", ".join(result_words))

find_words_with_non_first_syllable_stress()
<<<Шашлык, горелым, прибежал, сосед, дошёл, его>>>