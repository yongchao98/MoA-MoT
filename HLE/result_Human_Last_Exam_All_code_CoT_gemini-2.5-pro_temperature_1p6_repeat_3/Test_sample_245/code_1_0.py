import re

def find_words_with_non_first_syllable_accent():
    """
    Finds and prints words from a specific Russian text that have more than one syllable
    and are not accented on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary with lowercase words and their 0-indexed stressed syllable number.
    # This data is pre-analyzed as Russian stress is a linguistic feature.
    # A value of 0 means the first syllable is stressed.
    accent_data = {
        'шашлык': 1,
        'запах': 0,
        'горелым': 1,
        'вскоре': 0,
        'прибежал': 2,
        'сосед': 1,
        'видно': 0,
        'дошёл': 1, # The letter 'ё' is always stressed
        'его': 1,
        'дома': 0
    }
    
    vowels = "аеёиоуыэюя"

    # 1. Tokenize the text into words, preserving original case.
    # This regex finds all sequences of Cyrillic letters.
    words = re.findall(r'[а-яА-ЯёЁ]+', text)
    
    result_words = []
    
    for word in words:
        word_lower = word.lower()
        
        # 2. Filter out one-syllable words by counting vowels.
        num_syllables = sum(1 for char in word_lower if char in vowels)
        if num_syllables <= 1:
            continue
            
        # 3. Check if the accent is not on the first syllable (i.e., stress index > 0).
        # We look up the word in our pre-analyzed accent data.
        if word_lower in accent_data and accent_data[word_lower] > 0:
            result_words.append(word)

    # 4. Print the final list, comma-separated.
    print(", ".join(result_words))

find_words_with_non_first_syllable_accent()