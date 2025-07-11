import re

def find_words_with_non_initial_accent():
    """
    This function analyzes a given Russian sentence to find words that
    are multi-syllabic and do not have the accent on the first syllable.
    """
    # The source text in Russian
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # In Russian, a syllable is formed around a vowel.
    # "Accent on the first syllable" means the first vowel is stressed.
    # To solve this without a complex linguistics library, we'll store the accent
    # information for the relevant words.
    # The value is the 0-indexed position of the stressed vowel.
    # An accent on the first syllable corresponds to an index of 0.
    accent_data = {
        'шашлык': 1,   # ша-шлЫк (2nd vowel)
        'запах': 0,    # зА-пах (1st vowel)
        'горелым': 1,  # го-рЕ-лым (2nd vowel)
        'вскоре': 0,   # вскО-ре (1st vowel)
        'прибежал': 2, # при-бе-жАл (3rd vowel)
        'сосед': 1,    # со-сЕд (2nd vowel)
        'видно': 0,    # вИ-дно (1st vowel)
        'дошёл': 1,    # до-шЁл (2nd vowel, 'ё' is always stressed)
        'его': 1,      # е-гО (2nd vowel)
        'дома': 0      # дО-ма (genitive case from context, 1st vowel)
    }

    vowels = "аеёиоуыэюя"

    # Use regex to find all sequences of Cyrillic letters (words)
    words_in_text = re.findall(r'[а-яА-ЯёЁ]+', text)

    result_words = []
    # Iterate through the words in their order of appearance
    for word in words_in_text:
        word_lower = word.lower()

        # Count syllables by counting the number of vowels
        syllable_count = sum(1 for char in word_lower if char in vowels)

        # Condition 1: Exclude one-syllable words
        if syllable_count <= 1:
            continue

        # Condition 2: Check if the word's accent is not on the first syllable
        # (i.e., the stressed vowel index > 0)
        if word_lower in accent_data and accent_data[word_lower] > 0:
            result_words.append(word)

    # Print the final list, comma-separated.
    print(",".join(result_words))

find_words_with_non_initial_accent()
<<<Шашлык,горелым,прибежал,сосед,дошёл,его>>>