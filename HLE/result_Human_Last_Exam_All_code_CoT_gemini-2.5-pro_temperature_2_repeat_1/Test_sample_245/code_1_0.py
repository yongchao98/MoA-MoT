import re

def find_words_with_non_first_syllable_accent():
    """
    Analyzes a given Russian sentence to find words that are not one-syllable
    and do not have the accent on the first syllable.

    The accentuation for each word is pre-defined based on standard Russian.
    The combining acute accent mark is represented by `\\u0301`.
    """
    # The original Russian sentence from the prompt.
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A pre-defined mapping of words from the text to their correctly accented versions.
    # The `ё` is always stressed. Other stresses are marked with `\u0301`.
    # This dictionary helps look up the accent for each word found in the text.
    accent_map = {
        "шашлык": "шашлы\u0301к",
        "запах": "за\u0301пах",
        "горелым": "горе\u0301лым",
        "вскоре": "вско\u0301ре",
        "прибежал": "прибежа\u0301л",
        "сосед": "сосе\u0301д",
        "видно": "ви\u0301дно",
        "дошёл": "дошёл",
        "его": "его\u0301",
        "дома": "до\u0301ма",
    }
    
    # Russian vowels, used for counting syllables.
    vowels = "аеёиоуыэюя"
    
    # Split the sentence into potential words using spaces and common punctuation.
    words_in_order = re.split(r'[\s,\-—]+', text)
    
    result_words = []
    
    for word_orig in words_in_order:
        # Clean up any trailing punctuation from the word.
        word_clean = word_orig.strip(".,-— ")
        
        if not word_clean:
            continue

        # Convert to lowercase to match the keys in our accent map.
        word_lower = word_clean.lower()
        
        # Count syllables by counting the number of vowels.
        syllable_count = sum(1 for char in word_lower if char in vowels)
        
        # Rule: Exclude one-syllable words.
        if syllable_count <= 1:
            continue

        # Get the accented version of the word from our map.
        accented_word = accent_map.get(word_lower)
        if not accented_word:
            continue

        # Find the index of the first vowel.
        first_vowel_pos = -1
        for i, char in enumerate(accented_word):
            if char in vowels:
                first_vowel_pos = i
                break
        
        # Determine if the accent is on the first syllable.
        accent_is_on_first_syllable = False
        
        # Check for accent mark '\u0301' after the first vowel.
        if accented_word.find('\u0301') == first_vowel_pos + 1:
            accent_is_on_first_syllable = True
        # Check if the first vowel itself is 'ё'.
        elif accented_word[first_vowel_pos] == 'ё':
            accent_is_on_first_syllable = True
            
        # If the accent is not on the first syllable, add the original word to the list.
        if not accent_is_on_first_syllable:
            result_words.append(word_clean)
            
    # Print the final list, comma-separated.
    print(",".join(result_words))

find_words_with_non_first_syllable_accent()