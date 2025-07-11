import re

def find_words_with_non_first_syllable_accent():
    """
    This function finds and prints words from a specific Russian sentence
    that have more than one syllable and do not have the accent on the first syllable.
    """
    # The original text provided by the user.
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Based on phonetic analysis, this set contains the lowercase versions of the target words.
    # These are the multi-syllable words from the text where the accent is NOT on the first syllable.
    # Words excluded (accent on 1st syllable): запах, вскоре, видно, дома.
    # Words excluded (one syllable): и, до.
    target_words = {"шашлык", "горелым", "прибежал", "сосед", "дошёл", "его"}

    # Split the sentence into parts by spaces to maintain the order of words.
    words_in_order = text.split()

    # List to store the final results.
    result_list = []

    # Iterate through the words in their original order.
    for word in words_in_order:
        # Clean the word of any trailing punctuation (like commas, periods, hyphens).
        cleaned_word = re.sub(r'[^\w\s]', '', word)
        
        # Check if the cleaned, lowercase word is in our target set.
        if cleaned_word.lower() in target_words:
            result_list.append(cleaned_word)

    # Print the final list, with words separated by a comma.
    print(",".join(result_list))

find_words_with_non_first_syllable_accent()