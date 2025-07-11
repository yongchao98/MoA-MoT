import re

def find_words_with_non_first_syllable_stress():
    """
    Finds words from a given Russian text that are not monosyllabic and do not have the stress on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    
    # Pre-defined stress data for words in the text. The value represents the 1-indexed syllable that is stressed.
    # Note: 'ё' is always stressed, so its position determines the stress.
    stress_data = {
        "шашлык": 2,   # шаш-лЫк
        "запах": 1,    # зА-пах
        "горелым": 2,  # го-рЕ-лым
        "вскоре": 1,   # вскО-ре
        "прибежал": 3, # при-бе-жАл
        "сосед": 2,     # со-сЕд
        "видно": 1,     # вИд-но
        "дошёл": 2,      # до-шЁл
        "его": 2,      # е-гО
        "дома": 1,      # дО-ма
    }
    
    vowels = "аеёиоуыэюя"
    
    # Extract all words from the text, preserving order and case.
    words_in_order = re.findall(r'\b[а-яА-ЯёЁ]+\b', text)
    
    result_words = []
    
    for word in words_in_order:
        lower_word = word.lower()
        
        # Count syllables by counting vowels.
        syllable_count = sum(1 for char in lower_word if char in vowels)
        
        # Exclude one-syllable words.
        if syllable_count <= 1:
            continue
        
        # Check if the stress is not on the first syllable.
        # .get() is used to avoid an error if a word is not in our dictionary.
        stress_position = stress_data.get(lower_word)
        
        if stress_position and stress_position > 1:
            result_words.append(word)
            
    # Print the final comma-separated list.
    print(','.join(result_words))

find_words_with_non_first_syllable_stress()
<<<Шашлык,горелым,прибежал,сосед,дошёл,его>>>