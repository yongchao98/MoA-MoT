import re

def find_words_with_non_first_syllable_accent():
    """
    Finds and lists words from a given Russian text that are not monosyllabic
    and do not have the stress on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary mapping Russian words to their properties:
    # (number_of_syllables, index_of_stressed_syllable_starting_from_0)
    # This linguistic data is pre-analyzed.
    word_stress_data = {
        # word: (syllable_count, stressed_syllable_index)
        "шашлык": (2, 1),
        "запах": (2, 0),
        "горелым": (3, 1),
        "и": (1, 0),
        "вскоре": (2, 0),
        "прибежал": (3, 2),
        "сосед": (2, 1),
        "видно": (2, 0),
        "дошёл": (2, 1), # The letter "ё" is always stressed
        "до": (1, 0),
        "его": (2, 1),
        "дома": (2, 0)
    }

    # Extract all Cyrillic words from the text in order of appearance
    # using regular expressions.
    words_in_order = re.findall(r'\b[а-яА-ЯёЁ]+\b', text)

    result_list = []
    for word in words_in_order:
        # Use the lowercase version for dictionary lookup
        word_lower = word.lower()
        
        if word_lower in word_stress_data:
            syllables, stressed_syllable_index = word_stress_data[word_lower]
            
            # Check for two conditions:
            # 1. The word has more than one syllable.
            # 2. The stress is not on the first syllable (index > 0).
            if syllables > 1 and stressed_syllable_index > 0:
                result_list.append(word)

    # Join the collected words with a comma and a space for the final output.
    print(", ".join(result_list))

find_words_with_non_first_syllable_accent()