import re

def find_unstressed_words():
    """
    This function processes a given Russian sentence to find all multi-syllable words
    where the stress is not on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A pre-compiled dictionary mapping lowercase words to their stressed vowel's 0-index.
    # This is necessary as programmatic stress detection is a complex linguistic task.
    # For example, in 'шашлы́к', the stressed vowel 'ы' is the 2nd vowel, so its index is 1.
    # The letter 'ё' is always stressed.
    stress_data = {
        "шашлык": 1,    # шаш-лЫ́к
        "запах": 0,     # зА́-пах
        "горелым": 1,   # го-рЕ́-лым
        "вскоре": 0,    # вскО́-ре
        "прибежал": 2,  # при-бе-жА́л
        "сосед": 1,      # со-сЕ́д
        "видно": 0,      # вИ́д-но
        "дошёл": 1,     # до-шЁл
        "его": 1,       # е-гО́
        "дома": 0,      # дО́-ма (in this context)
    }

    # Russian vowels used for syllable counting.
    vowels = "аеёиоуыэюя"

    # Use regex to extract all words, preserving their original form.
    original_words = re.findall(r'[\w-]+', text)
    
    result_list = []
    
    for word in original_words:
        word_lower = word.lower()
        
        # Count syllables by counting vowels.
        syllable_count = sum(1 for char in word_lower if char in vowels)
        
        # Process only multi-syllable words present in our stress data.
        if syllable_count > 1 and word_lower in stress_data:
            # Check if the stress is NOT on the first syllable (index > 0).
            stressed_vowel_index = stress_data[word_lower]
            if stressed_vowel_index > 0:
                result_list.append(word)

    # Print the final list, comma-separated.
    print(",".join(result_list))

find_unstressed_words()
<<<Шашлык,горелым,прибежал,сосед,дошёл,его>>>