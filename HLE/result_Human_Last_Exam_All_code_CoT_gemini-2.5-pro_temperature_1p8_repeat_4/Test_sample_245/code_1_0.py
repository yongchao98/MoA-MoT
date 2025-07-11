import re

def solve():
    """
    Finds and prints Russian words from a given text that have more than one
    syllable and are not stressed on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary mapping lowercase multi-syllable words to their stressed vowel's index (0-based).
    # For example, in "горелым" (го-РЕ-лым), the vowels are 'о', 'е', 'ы'. Stress is on 'е' (index 1).
    # Words with stress on the first syllable (e.g., "запах" - ЗА-пах) have an index of 0.
    stress_data = {
        'шашлык': 1,
        'запах': 0,
        'горелым': 1,
        'вскоре': 0,
        'прибежал': 2,
        'сосед': 1,
        'видно': 0,
        'дошёл': 1,  # Stress on 'ё' is always the case.
        'его': 1,
        'дома': 0
    }

    vowels = "аеёиоуыэюя"
    
    # Use a regular expression to find all words (sequences of Cyrillic letters) in order.
    words_in_order = re.findall(r'[а-яА-ЯёЁ]+', text)
    
    result_words = []
    
    for word in words_in_order:
        word_lower = word.lower()
        
        # Count vowels to determine the number of syllables.
        syllable_count = sum(1 for char in word_lower if char in vowels)
        
        # Rule 1: Exclude one-syllable words.
        if syllable_count <= 1:
            continue
            
        # Rule 2: Check if the word's stress is NOT on the first syllable.
        # We check if the word is in our data and its stress index is greater than 0.
        stressed_vowel_index = stress_data.get(word_lower)
        if stressed_vowel_index is not None and stressed_vowel_index > 0:
            result_words.append(word)

    # Print the final list, comma-separated.
    print(", ".join(result_words))

solve()