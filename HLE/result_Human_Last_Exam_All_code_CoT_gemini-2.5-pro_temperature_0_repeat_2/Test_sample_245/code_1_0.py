import re

def solve():
    """
    Finds and prints words from a Russian text that are not monosyllabic
    and do not have the accent on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary mapping words from the text to their accent position.
    # The accent position is 1-indexed (1 for the first syllable, 2 for the second, etc.).
    # This is necessary because Russian accent is lexical and not rule-based.
    # Words with accent on the first syllable: запах, видно, дома.
    # Words with accent not on the first syllable: шашлык, горелым, вскоре, прибежал, сосед, дошёл, его.
    accented_words = {
        "шашлык", "горелым", "вскоре", "прибежал", "сосед", "дошёл", "его"
    }

    # A set of one-syllable words to exclude.
    one_syllable_words = {"и", "до"}

    # Use regex to extract all Russian words, preserving order.
    words = re.findall(r'\b[а-яА-ЯёЁ]+\b', text)

    result_list = []
    for word in words:
        # Normalize the word to lowercase for checking.
        lower_word = word.lower()

        # Check if the word is a known multi-syllable word with the accent
        # not on the first syllable.
        if lower_word not in one_syllable_words and lower_word in accented_words:
            result_list.append(word)

    # Print the final list, comma-separated.
    print(", ".join(result_list))

solve()