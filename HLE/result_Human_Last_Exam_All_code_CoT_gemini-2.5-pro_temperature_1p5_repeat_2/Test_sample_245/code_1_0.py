import re

def find_words_with_non_initial_accent():
    """
    Finds and prints words from a given Russian text that are not one-syllable 
    and do not have the accent on the first syllable.
    """
    # The original Russian text
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Based on linguistic analysis, these are the target words in lowercase.
    # They have more than one syllable and the stress is not on the first syllable.
    # Шашлы́к, горе́лым, прибежа́л, сосе́д, дошёл, его́
    target_words_set = {'шашлык', 'горелым', 'прибежал', 'сосед', 'дошёл', 'его'}

    # Clean punctuation and split the text into words, preserving original casing.
    # We replace punctuation with spaces to ensure correct splitting.
    cleaned_text = re.sub(r'[,.-]', ' ', text)
    all_words = cleaned_text.split()
    
    # Collect the words that match our criteria in their order of appearance.
    result_words = []
    for word in all_words:
        if word.lower() in target_words_set:
            result_words.append(word)
            
    # Print the final list, comma-separated.
    print(", ".join(result_words))

find_words_with_non_initial_accent()