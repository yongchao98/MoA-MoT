import re

def find_words_with_non_first_syllable_stress():
    """
    Finds and prints words from a given Russian text that are not single-syllable
    and do not have the stress on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary mapping lowercase words to their 0-indexed stressed syllable.
    # This information is based on linguistic knowledge of Russian.
    stress_data = {
        "шашлык": 1,   # шаш-лы́к
        "запах": 0,    # за́-пах
        "горелым": 1,  # го-ре́-лым
        "вскоре": 0,   # вско́-ре
        "прибежал": 2, # при-бе-жа́л
        "сосед": 1,    # со-се́д
        "видно": 0,    # ви́д-но
        "дошёл": 1,    # до-шёл
        "его": 1,      # е-го́
        "дома": 0,     # до́-ма
    }

    vowels = "аеёиоуыэюя"
    
    # Extract all sequences of Cyrillic letters (words) from the text.
    words = re.findall(r'[а-яА-ЯёЁ]+', text)

    result_words = []
    
    for word in words:
        word_lower = word.lower()
        
        # Count the number of syllables by counting vowels.
        syllable_count = sum(1 for char in word_lower if char in vowels)
        
        # Rule: Exclude one-syllable words.
        if syllable_count <= 1:
            continue
        
        stressed_syllable_index = -1
        
        # Rule: The letter 'ё' is always stressed.
        if 'ё' in word_lower:
            vowel_counter = 0
            for char in word_lower:
                if char in vowels:
                    if char == 'ё':
                        stressed_syllable_index = vowel_counter
                        break
                    vowel_counter += 1
        # Otherwise, look up the stress position in our data.
        elif word_lower in stress_data:
            stressed_syllable_index = stress_data[word_lower]
            
        # Add the word to our result list if stress is not on the first syllable (index > 0).
        if stressed_syllable_index > 0:
            result_words.append(word)
            
    # Print the final list, comma-separated.
    print(", ".join(result_words))

find_words_with_non_first_syllable_stress()
<<<Шашлык, горелым, прибежал, сосед, дошёл, его>>>