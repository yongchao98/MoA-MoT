import re

def find_words_with_non_initial_stress():
    """
    Analyzes a Russian sentence to find and list all multi-syllable words
    where the stress does not fall on the first syllable.
    """
    # The source text provided by the user.
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary is used to provide stress information for words from the text.
    # The key is the lowercase word, and the value is the word with an accent
    # mark (´) after the stressed vowel, or containing 'ё', which is always stressed.
    stress_data = {
        "шашлык": "шашлы́к",
        "запах": "за́пах",
        "горелым": "горе́лым",
        "вскоре": "вско́ре",
        "прибежал": "прибежа́л",
        "сосед": "сосе́д",
        "видно": "ви́дно",
        "дошёл": "дошёл",
        "его": "его́",
        "дома": "до́ма"
    }

    # A string containing all Russian vowels for syllable counting.
    vowels = "аеёиоуыэюя"

    # Use regular expressions to extract all word-like sequences from the text,
    # converting them to lowercase for consistent matching.
    words_in_order = re.findall(r'[\w-]+', text.lower())

    result_words = []
    for word in words_in_order:
        # 1. Filter out words with one syllable (i.e., one vowel).
        syllable_count = sum(1 for char in word if char in vowels)
        if syllable_count <= 1:
            continue

        # Proceed only if we have stress data for the word.
        if word in stress_data:
            stressed_word_representation = stress_data[word]

            # 2. Determine if the stress is on the first syllable.
            # First, find the index of the first vowel in the word.
            first_vowel_index = -1
            for i, char in enumerate(word):
                if char in vowels:
                    first_vowel_index = i
                    break

            # A word has stress on the first syllable if its first vowel is 'ё'
            # or if it's marked with an accent in our stress data.
            is_stress_on_first = False
            if word[first_vowel_index] == 'ё':
                is_stress_on_first = True
            # The combining accent mark '´' appears *after* the character it modifies.
            elif (first_vowel_index + 1 < len(stressed_word_representation) and
                  stressed_word_representation[first_vowel_index + 1] == '́'):
                is_stress_on_first = True

            # 3. If the stress is NOT on the first syllable, add it to the results.
            if not is_stress_on_first:
                result_words.append(word)

    # 4. Print the final list as a comma-separated string.
    print(",".join(result_words))

# Run the analysis and print the result.
find_words_with_non_initial_stress()