import re

def find_words_with_non_initial_stress():
    """
    Finds and prints Russian words from a given text that have more than one syllable
    and the stress is not on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A pre-defined dictionary with stress information for multi-syllable words in the text.
    # The value is the 0-indexed position of the stressed vowel.
    # Шашлы́к (a,ы́ -> 1), за́пах (а́,а -> 0), горе́лым (о,е́,ы -> 1), вско́ре (о́,е -> 0),
    # прибежа́л (и,е,а́ -> 2), сосе́д (о,е́ -> 1), ви́дно (и́,о -> 0), дошё́л (о,ё́ -> 1),
    # его́ (е,о́ -> 1), до́ма (о́,а -> 0).
    stress_data = {
        'шашлык': 1,
        'запах': 0,
        'горелым': 1,
        'вскоре': 0,
        'прибежал': 2,
        'сосед': 1,
        'видно': 0,
        'дошёл': 1,
        'его': 1,
        'дома': 0
    }

    vowels = "аоуыэяеёюи"
    
    # Use regex to find all sequences of Russian letters (words)
    words = re.findall(r'[а-яА-ЯёЁ]+', text)
    
    result_words = []
    
    for word in words:
        # Normalize to lowercase for dictionary lookup and vowel counting
        lower_word = word.lower()
        
        # Count vowels to approximate syllables
        vowel_count = sum(1 for char in lower_word if char in vowels)
        
        # Condition 1: Must be a multi-syllable word
        if vowel_count > 1:
            # Check if we have stress information for this word
            if lower_word in stress_data:
                # Condition 2: Stress is not on the first syllable (stressed vowel index > 0)
                if stress_data[lower_word] > 0:
                    result_words.append(word)

    print(",".join(result_words))

find_words_with_non_initial_stress()