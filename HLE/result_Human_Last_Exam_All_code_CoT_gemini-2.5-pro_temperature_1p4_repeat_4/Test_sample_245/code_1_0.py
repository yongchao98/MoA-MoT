import re

def find_words_with_non_first_syllable_stress():
    """
    Analyzes a Russian text to find words with more than one syllable
    where the stress is not on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    
    # Russian vowels used for syllable counting.
    vowels = "аеёиоуыэюя"
    
    # A dictionary mapping each lowercase multi-syllable word to the 0-based index of its stressed vowel.
    # For example, in "горелым", the vowels are (о, е, ы) and the stress is on 'е', so the index is 1.
    # The 'ё' is always stressed, so for "дошёл" (vowels о, ё), the stress index is 1.
    stress_data = {
        "шашлык": 1,
        "запах": 0,
        "горелым": 1,
        "вскоре": 0,
        "прибежал": 2,
        "сосед": 1,
        "видно": 0,
        "дошёл": 1,
        "его": 1,
        "дома": 0
    }
    
    # 1. Tokenize the text into words, preserving original capitalization.
    original_words = re.findall(r'[а-яА-ЯёЁ]+', text)
    
    result_words = []
    
    # 3. Filter the words based on the criteria.
    for word in original_words:
        # Use the lowercase version for analysis.
        word_lower = word.lower()
        
        # a. Count syllables by counting vowels.
        num_vowels = sum(1 for char in word_lower if char in vowels)
        
        # Skip one-syllable words.
        if num_vowels <= 1:
            continue
            
        # b. Look up the stress position.
        # We assume all necessary words are in our stress_data dictionary.
        if word_lower in stress_data:
            stressed_vowel_index = stress_data[word_lower]
            
            # c. & d. If the stress is not on the first syllable (index > 0), keep the word.
            if stressed_vowel_index > 0:
                result_words.append(word)

    # 4. Format and print the final output.
    print(",".join(result_words))

find_words_with_non_first_syllable_stress()
<<<Шашлык,горелым,прибежал,сосед,дошёл,его>>>