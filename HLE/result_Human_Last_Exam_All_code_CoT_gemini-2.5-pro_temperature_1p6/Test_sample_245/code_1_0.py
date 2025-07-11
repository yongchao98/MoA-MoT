import re

def find_words_with_non_first_syllable_accent():
    """
    This function processes a Russian sentence to find and list all words
    (excluding one-syllable words) that do not have the accent on the first syllable.

    The logic relies on a pre-identified list of target words, as programmatic
    accent determination is complex. The script finds these words in the original
    sentence to maintain the order of appearance.
    """
    # The source text in Russian
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Manually identified list of target words that fit the criteria:
    # multi-syllable and accent not on the first syllable.
    # Stored in a set for efficient lookup.
    # шашлы́к, горе́лым, прибежа́л, сосе́д, дошёл, его́
    target_words = {"шашлык", "горелым", "прибежал", "сосед", "дошёл", "его"}

    # Use regex to find all word-like sequences from the lowercased text.
    # This correctly handles splitting the sentence into words.
    words = re.findall(r'\b[\w-]+\b', text.lower())

    # Collect the words from the sentence that are in our target list,
    # preserving their original order of appearance.
    result_list = []
    for word in words:
        if word in target_words:
            result_list.append(word)

    # Print the final list, with words separated by commas.
    print(",".join(result_list))

find_words_with_non_first_syllable_accent()