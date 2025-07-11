import re

def find_unstressed_first_syllable_words():
    """
    Finds and prints words from a specific Russian text that have more than one syllable
    and where the stress is not on the first syllable.
    """
    # The input text in Russian
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # This dictionary provides the necessary linguistic information about word stress.
    # The key is the word, and the value is the 0-based index of the stressed vowel.
    # An index of 0 means the first syllable is stressed.
    # An index > 0 means a later syllable is stressed.
    stress_data = {
        'шашлык': 1,   # ша-шлы́к (2nd vowel is stressed)
        'запах': 0,    # за́-пах (1st vowel is stressed)
        'горелым': 1,  # го-ре́-лым (2nd vowel is stressed)
        'вскоре': 0,   # вско́-ре (1st vowel is stressed)
        'прибежал': 2, # при-бе-жа́л (3rd vowel is stressed)
        'сосед': 1,     # со-се́д (2nd vowel is stressed)
        'видно': 0,     # ви́д-но (1st vowel is stressed)
        'дошёл': 1,    # до-шёл (2nd vowel is stressed)
        'его': 1,      # е-го́ (2nd vowel is stressed)
        'дома': 0      # до́-ма (genitive singular of 'дом', 1st vowel is stressed)
    }

    # A string containing all Russian vowels for syllable counting
    vowels = "аеёиоуыэюя"

    # Use regex to find all sequences of Cyrillic letters, and convert to lowercase
    words = re.findall(r'[а-яё]+', text.lower())

    # A list to store the words that match our criteria
    result_list = []

    # Iterate through the words in the order they appeared in the text
    for word in words:
        # Count the number of syllables by counting the vowels
        syllable_count = sum(1 for char in word if char in vowels)

        # Skip one-syllable words
        if syllable_count <= 1:
            continue

        # Get the stressed vowel index from our data. Default to -1 if not found.
        stressed_vowel_index = stress_data.get(word, -1)

        # If the stressed vowel's index is greater than 0, the stress is not on the first syllable.
        if stressed_vowel_index > 0:
            result_list.append(word)

    # Join the results into a single comma-separated string and print it
    print(",".join(result_list))

find_unstressed_first_syllable_words()