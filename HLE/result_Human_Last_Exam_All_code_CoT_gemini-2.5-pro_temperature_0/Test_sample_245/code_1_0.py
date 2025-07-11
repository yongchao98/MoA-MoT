import re

def find_and_print_words():
    """
    This function identifies and prints words from a specific Russian sentence
    that are multi-syllabic and do not have the stress on the first syllable.
    """
    # The source text in Russian
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Based on linguistic analysis, these are the lowercase words that meet the criteria:
    # - More than one syllable.
    # - Stress is not on the first syllable.
    target_words = {"шашлык", "горелым", "прибежал", "сосед", "дошёл", "его"}

    # Use regex to extract all words (sequences of Cyrillic letters) from the text,
    # preserving their original order and case.
    words_in_order = re.findall(r'[а-яА-ЯёЁ]+', text)

    # Collect the words that match our criteria
    result_list = []
    for word in words_in_order:
        # Check if the lowercase version of the word is in our target set
        if word.lower() in target_words:
            result_list.append(word)

    # Print the final list, comma-separated, as requested.
    print(", ".join(result_list))

# Execute the function to print the result
find_and_print_words()