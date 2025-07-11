import re

def solve():
    """
    Finds and prints words from a given Russian text that are not monosyllabic
    and do not have the stress on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Manually determined stress data: {lowercase_word: (syllable_count, stressed_syllable_number)}
    # The number of syllables is based on the number of vowels.
    # The stress position is based on standard Russian pronunciation.
    stress_info = {
        'шашлык':   (2, 2),  # Шаш-лы́к
        'запах':    (2, 1),  # за́-пах
        'горелым':  (3, 2),  # го-ре́-лым
        'и':        (1, 1),  # и
        'вскоре':   (2, 1),  # вско́-ре
        'прибежал': (3, 3),  # при-бе-жа́л
        'сосед':     (2, 2),  # со-се́д
        'видно':     (2, 1),  # ви́д-но
        'дошёл':    (2, 2),  # до-шёл
        'до':       (1, 1),  # до
        'его':      (2, 2),  # е-го́
        'дома':     (2, 1)   # до́-ма
    }

    # Extract all words from the text using regex
    words_in_text = re.findall(r'\b[а-яА-ЯёЁ]+\b', text)

    result_list = []
    for word in words_in_text:
        # Look up the word's properties in our analysis table
        word_data = stress_info.get(word.lower())

        if word_data:
            syllable_count, stressed_syllable = word_data

            # Condition 1: Exclude one-syllable words
            is_multisyllabic = syllable_count > 1

            # Condition 2: Exclude words with stress on the first syllable
            is_stress_not_on_first = stressed_syllable != 1

            if is_multisyllabic and is_stress_not_on_first:
                result_list.append(word)

    print(",".join(result_list))

solve()
<<<Шашлык,горелым,прибежал,сосед,дошёл,его>>>