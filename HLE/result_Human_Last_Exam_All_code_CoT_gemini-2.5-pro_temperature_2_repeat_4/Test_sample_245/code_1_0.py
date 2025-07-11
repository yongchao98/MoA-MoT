import re

def find_words_with_non_initial_stress():
    """
    Finds and prints words from a Russian sentence that are multi-syllable
    and do not have stress on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary mapping a lowercase word to the character index of its stressed vowel.
    # This linguistic data is pre-determined.
    stress_data_by_index = {
        'шашлык': 4,
        'запах': 1,
        'горелым': 3,
        'вскоре': 3,
        'прибежал': 7,
        'сосед': 4,
        'видно': 1,
        'его': 2,
        'дома': 1
    }
    
    vowels = "аеёиоуыэюя"
    
    # Clean punctuation and split into a list of words.
    clean_text = re.sub(r'[^\w\s]', '', text)
    words = clean_text.split()
    
    result_words = []
    
    for word in words:
        lower_word = word.lower()
        
        # Find the indices of all vowels in the word.
        vowel_indices = []
        for i, char in enumerate(lower_word):
            if char in vowels:
                vowel_indices.append(i)
                
        # Skip one-syllable words.
        if len(vowel_indices) <= 1:
            continue
            
        stressed_vowel_index = -1
        # The letter 'ё' is always stressed.
        if 'ё' in lower_word:
            stressed_vowel_index = lower_word.find('ё')
        elif lower_word in stress_data_by_index:
            stressed_vowel_index = stress_data_by_index[lower_word]
        
        # Check if the word's stress is not on the first syllable.
        if stressed_vowel_index != -1:
            first_vowel_index = vowel_indices[0]
            if stressed_vowel_index != first_vowel_index:
                result_words.append(word)

    print(','.join(result_words))

find_words_with_non_initial_stress()