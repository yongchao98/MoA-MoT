import re

def find_words_with_non_first_syllable_accent():
    """
    Analyzes a given Russian text to find multi-syllable words where the
    accent is not on the first syllable.
    """
    # The text to analyze
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A pre-compiled dictionary of the multi-syllable words from the text.
    # The key is the lowercase word, and the value is the 1-based index of the stressed syllable.
    # One-syllable words (like 'и', 'до') are excluded from this dictionary.
    stress_data = {
        "шашлык": 2,    # шаш-лы́к
        "запах": 1,     # за́-пах
        "горелым": 2,   # го-ре́-лым
        "вскоре": 1,    # вско́-ре
        "прибежал": 3,  # при-бе-жа́л
        "сосед": 2,     # со-се́д
        "видно": 1,     # ви́д-но
        "дошёл": 2,     # до-шёл (the vowel 'ё' is always stressed)
        "его": 2,       # е-го́
        "дома": 1       # до́-ма
    }

    # Split the text into words using spaces and punctuation as delimiters
    # to maintain the original order.
    words = re.split(r'[\s,\-]+', text)

    # A list to store the words that match the criteria.
    result_words = []

    # Iterate through each word in its order of appearance.
    for word in words:
        # Clean the word by removing trailing punctuation and converting to lowercase for lookup.
        clean_word = word.strip('.,').lower()

        # Check if the word is a multi-syllable word we have stress data for.
        if clean_word in stress_data:
            # Check if the stress is not on the first syllable.
            if stress_data[clean_word] > 1:
                # If it matches, add the original word (cleaned of punctuation) to our results.
                result_words.append(word.strip('.,'))

    # Print the final list, comma-separated, as requested.
    print(", ".join(result_words))

# Execute the function to print the result.
find_words_with_non_first_syllable_accent()