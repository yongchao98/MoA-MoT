import re

def find_unstressed_words():
    """
    Finds and prints words from a given Russian text that are not stressed on the first syllable
    and are not one-syllable words.
    """
    # The input text provided by the user.
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Russian accent/stress is lexical. We need to know the stress for each word.
    # This dictionary maps a lowercase word to the 0-based index of its stressed vowel.
    # For example, в 'запах' (зА-пах), vowels are 'а' and 'а'. The stress is on the first 'а', so its index is 0.
    # In 'шашлык' (ша-шлЫк), vowels are 'а' and 'ы'. The stress is on 'ы', the second vowel, so its index is 1.
    stress_data = {
        "шашлык": 1,    # stress on ы (2nd vowel)
        "запах": 0,     # stress on first а (1st vowel)
        "горелым": 1,   # stress on е (2nd vowel)
        "вскоре": 0,    # stress on о (1st vowel)
        "прибежал": 2,  # stress on second а (3rd vowel)
        "сосед": 1,    # stress on е (2nd vowel)
        "видно": 0,     # stress on и (1st vowel)
        "дошёл": 1,     # stress on ё (2nd vowel)
        "его": 1,       # stress on о (2nd vowel)
        "дома": 0       # stress on о (1st vowel, from "до его дОма")
    }

    # The set of Russian vowel characters.
    vowels = "аеёиоуыэюя"

    # Clean the text by replacing punctuation with spaces, then split into words.
    # This preserves the original order and casing of the words.
    cleaned_text = re.sub(r'[,\-.]', ' ', text)
    words_in_order = cleaned_text.split()

    result_words = []
    for word in words_in_order:
        word_lower = word.lower()

        # Count vowels to determine the number of syllables.
        syllable_count = sum(1 for char in word_lower if char in vowels)

        # Condition 1: Exclude one-syllable words.
        if syllable_count <= 1:
            continue

        # Check if we have stress information for this word.
        if word_lower in stress_data:
            stressed_vowel_index = stress_data[word_lower]
            # Condition 2: Check that the stress is NOT on the first syllable (index 0).
            if stressed_vowel_index > 0:
                result_words.append(word)

    # Print the final list, comma-separated, as requested.
    print(",".join(result_words))

find_unstressed_words()