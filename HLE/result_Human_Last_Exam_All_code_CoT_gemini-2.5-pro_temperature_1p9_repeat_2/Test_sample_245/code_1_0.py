import re

def find_words_with_non_first_syllable_accent():
    """
    This function processes a Russian text to find all words that are not
    one-syllable and do not have the stress on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Since determining Russian word stress programmatically is complex and requires
    # a comprehensive dictionary, we will pre-analyze the words in the specific text.
    # The following set contains the lowercase versions of the target words that
    # are multi-syllable and have their accent on a syllable other than the first.
    target_words_lower = {
        "шашлык",   # шаш-лы́к
        "горелым",  # го-ре́-лым
        "прибежал", # при-бе-жа́л
        "сосед",     # со-се́д
        "дошёл",    # до-шёл (ё is always stressed)
        "его"       # е-го́
    }

    # Use regex to find all word-like sequences in the original text to maintain order.
    # This captures words containing hyphens as well.
    words_in_text = re.findall(r'[\w-]+', text)

    # Filter the list of words from the text, keeping only our targets.
    # The check is case-insensitive, but we keep the original capitalization.
    result_list = [word for word in words_in_text if word.lower() in target_words_lower]

    # Print the final list of words, separated by commas.
    print(",".join(result_list))

find_words_with_non_first_syllable_accent()