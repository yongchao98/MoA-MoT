# Before running, you may need to install the required library.
# You can do this by running the following command in your terminal:
# pip install russian-tc

import re
import russian_tc

def find_words_with_non_initial_stress(text):
    """
    Finds words in a Russian text that are not one-syllable and do not have stress on the first syllable.
    """
    # Using a set for fast lookups
    vowels = "аеёиоуыэюя"
    
    # Extract all words, converting to lowercase. [а-яё]+ matches sequences of Russian letters.
    # We use re.findall to maintain the original order of words.
    words = re.findall(r'[а-яё]+', text.lower())
    
    result_words = []

    for word in words:
        # Count vowels to identify polysyllabic words
        vowel_count = sum(1 for char in word if char in vowels)
        if vowel_count <= 1:
            continue

        first_vowel_index = -1
        for i, char in enumerate(word):
            if char in vowels:
                first_vowel_index = i
                break
        
        stressed_vowel_index = -1
        
        # The letter 'ё' is always stressed
        if 'ё' in word:
            stressed_vowel_index = word.find('ё')
        else:
            # Use the library to find the stress
            accented_word = russian_tc.accent(word)
            if "'" in accented_word:
                # The accent mark (') is placed after the stressed vowel.
                # Its index in the accented word corresponds to the index of the stressed vowel
                # in the original word.
                stressed_vowel_index = accented_word.find("'") - 1
        
        # We only add the word if the stress is found and is not on the first syllable.
        # The stressed vowel's index must be greater than the first vowel's index.
        if stressed_vowel_index > first_vowel_index:
            # We use the original word from our list, not the accented one
            result_words.append(word)
            
    # Print the final result as a comma-separated string
    print(", ".join(result_words))

# The input text from the user
text_from_user = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

# Run the function with the user's text
find_words_with_non_initial_stress(text_from_user)
<<<шашлык, горелым, прибежал, сосед, дошёл, его>>>