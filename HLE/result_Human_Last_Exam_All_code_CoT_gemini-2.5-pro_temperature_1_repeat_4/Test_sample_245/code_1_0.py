import re

def find_unstressed_words():
    """
    This function analyzes a given Russian sentence to find multi-syllable words
    where the stress does not fall on the first syllable.
    """
    # The sentence in Russian
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Russian vowels used for syllable counting
    vowels = "аеёиоуыэюя"

    # A pre-defined dictionary mapping each relevant word (lowercase) to the
    # zero-based index of its stressed vowel.
    # For example, in "шашлык" (шаш-ЛЫК), the vowels are 'а' and 'ы'.
    # The stress is on 'ы', which is the second vowel, so its index is 1.
    stress_data = {
        "шашлык": 1,   # шашлЫк
        "запах": 0,    # зАпах
        "горелым": 1,  # горЕлым
        "вскоре": 0,   # вскОре
        "прибежал": 2, # прибежАл
        "сосед": 1,    # сосЕд
        "видно": 0,    # вИдно
        "дошёл": 1,    # дошЁл ('ё' is the second vowel, so index 1)
        "его": 1,      # егО
        "дома": 0      # дОма
    }

    # Use regex to find all sequences of Cyrillic letters, preserving original case.
    # This extracts words and ignores punctuation.
    words = re.findall(r'\b[а-яА-ЯёЁ]+\b', text)

    result_words = []
    for word in words:
        # Use a lowercase version for analysis and dictionary lookup
        clean_word = word.lower()

        # Count syllables by counting the number of vowels
        syllable_count = sum(1 for char in clean_word if char in vowels)

        # We are only interested in words with more than one syllable
        if syllable_count > 1:
            # Check if the word is in our stress dictionary
            if clean_word in stress_data:
                # Get the index of the stressed vowel (0 for 1st, 1 for 2nd, etc.)
                stressed_vowel_index = stress_data[clean_word]

                # If the stressed vowel index is greater than 0, the stress is NOT
                # on the first syllable.
                if stressed_vowel_index > 0:
                    result_words.append(word)

    # Print the final list of words, comma-separated, as requested.
    print(",".join(result_words))

find_unstressed_words()