import re

def find_unstressed_words():
    """
    This function analyzes a Russian text to find words that are not
    one-syllable and do not have the accent on the first syllable.
    """
    # The text to analyze
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary mapping lowercase words to the 0-based index of their stressed vowel.
    # This is a practical approach for this specific problem, as determining
    # Russian stress programmatically is complex.
    stress_data = {
        'шашлык': 1,  # шашлЫк (vowels а, ы -> 2nd vowel)
        'запах': 0,   # зАпах (vowels а, а -> 1st vowel)
        'горелым': 1, # горЕлым (vowels о, е, ы -> 2nd vowel)
        'вскоре': 0,  # вскОре (vowels о, е -> 1st vowel)
        'прибежал': 2,# прибежАл (vowels и, е, а -> 3rd vowel)
        'сосед': 1,  # сосЕд (vowels о, е -> 2nd vowel)
        'видно': 0,   # вИдно (vowels и, о -> 1st vowel)
        'его': 1,     # егО (vowels е, о -> 2nd vowel)
        'дома': 0     # дОма (genitive case) (vowels о, а -> 1st vowel)
    }

    # Set of Russian vowels
    vowels = "аеёиоуыэюя"

    # Clean punctuation (commas, hyphens, periods) and split the text into words
    words = re.sub(r'[,\-\.]', '', text).split()

    result_words = []

    for word in words:
        # Create a list of vowels from the lowercase version of the word
        word_vowels = [char for char in word.lower() if char in vowels]

        # Rule 1: Exclude one-syllable words (words with 1 or 0 vowels)
        if len(word_vowels) <= 1:
            continue

        stressed_vowel_index = -1
        word_lower = word.lower()
        
        # The letter 'ё' is always stressed.
        if 'ё' in word_lower:
            stressed_vowel_index = word_vowels.index('ё')
        # Otherwise, look up the stress position in our predefined data.
        elif word_lower in stress_data:
            stressed_vowel_index = stress_data[word_lower]

        # Rule 2: The accent must NOT be on the first syllable (index 0).
        # We check for index > 0.
        if stressed_vowel_index > 0:
            result_words.append(word)

    # Print the final result, comma-separated, in the order of appearance.
    print(','.join(result_words))

find_unstressed_words()