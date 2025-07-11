import re

def solve_russian_stress():
    """
    Finds and prints multi-syllable Russian words from a given text
    that are not stressed on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # In Russian linguistics, stress determination is complex. For this specific
    # text, we will use a pre-defined dictionary mapping each relevant word to its
    # stressed syllable's index (0-based).
    # The letter 'ё' is always stressed, so we can identify its position.
    # Word: stressed syllable index
    stress_data = {
        'шашлык': 1,   # шаш-лЫк
        'запах': 0,    # зА-пах
        'горелым': 1,  # го-рЕ-лым
        'вскоре': 0,   # вскО-ре
        'прибежал': 2, # при-бе-жАл
        'сосед': 1,   # со-сЕд
        'видно': 0,   # вИд-но
        'дошёл': 1,   # до-шЁл (special case with 'ё')
        'его': 1,      # е-гО
        'дома': 0      # дО-ма
    }

    vowels = "аеиоуыэюяё"
    
    # Use regex to find all sequences of Cyrillic letters
    words = re.findall(r'[а-яёА-ЯЁ]+', text)
    
    result_words = []
    
    for word in words:
        word_lower = word.lower()
        
        # Count syllables by counting vowels
        syllable_count = sum(1 for char in word_lower if char in vowels)
        
        # We are only interested in words with more than one syllable
        if syllable_count > 1:
            # Check if the stress is NOT on the first syllable using our pre-defined data
            if word_lower in stress_data and stress_data[word_lower] > 0:
                result_words.append(word_lower)
    
    print(",".join(result_words))

solve_russian_stress()
<<<шашлык,горелым,прибежал,сосед,дошёл,его>>>