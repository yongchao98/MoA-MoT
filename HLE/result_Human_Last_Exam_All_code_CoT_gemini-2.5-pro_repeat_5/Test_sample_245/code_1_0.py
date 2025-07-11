import re

def find_words_with_non_initial_stress():
    """
    Analyzes a Russian sentence to find multi-syllable words where the stress
    is not on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A pre-analyzed dictionary for words in the text.
    # Key: lowercase word
    # Value: a tuple (total_syllables, stressed_syllable_number)
    stress_info = {
        "шашлык": (2, 2),
        "запах": (2, 1),
        "горелым": (3, 2),
        "вскоре": (2, 1),
        "прибежал": (3, 3),
        "сосед": (2, 2),
        "видно": (2, 1),
        "дошёл": (2, 2),
        "его": (2, 2),
        "дома": (2, 1)
    }

    # Extract all Russian words from the text, maintaining order.
    words = re.findall(r'[а-яА-ЯёЁ]+', text)

    result_words = []
    for word in words:
        # Use the lowercase version for dictionary lookup.
        word_lower = word.lower()

        # Check if we have stress information for this word.
        if word_lower in stress_info:
            syllable_count, stressed_syllable = stress_info[word_lower]

            # Apply the conditions: more than 1 syllable AND stress not on the 1st.
            if syllable_count > 1 and stressed_syllable > 1:
                # Add the original word to our results list.
                result_words.append(word)

    # Print the final list, comma-separated.
    print(",".join(result_words))

find_words_with_non_initial_stress()