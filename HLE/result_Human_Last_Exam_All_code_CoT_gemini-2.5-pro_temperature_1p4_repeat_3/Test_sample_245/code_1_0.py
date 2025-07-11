import re

def solve_russian_stress_task():
    """
    Finds words in a Russian sentence that are not one-syllable and
    do not have the stress on the first syllable.
    """
    # The input text
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Step 1: Pre-analyzed linguistic data.
    # The key is the lowercase word, the value is a tuple: (number_of_syllables, stressed_syllable_index)
    # The stressed syllable index is 1-based.
    word_data = {
        "шашлык": (2, 2),   # шаш-ЛЫК
        "запах": (2, 1),    # ЗА-пах
        "горелым": (3, 2),  # го-РЕ-лым
        "и": (1, 1),
        "вскоре": (2, 1),   # ВСКО-ре
        "прибежал": (3, 3), # при-бе-ЖАЛ
        "сосед": (2, 2),    # со-СЕД
        "видно": (2, 1),    # ВИД-но
        "дошёл": (2, 2),    # до-ШЁЛ (the letter ё is always stressed)
        "до": (1, 1),
        "его": (2, 2),      # е-ГО
        "дома": (2, 1),     # ДО-ма
    }

    # Step 2: Tokenize the sentence and remove punctuation
    # We use regex to find all sequences of Cyrillic letters.
    words = re.findall(r'\b[а-яА-ЯёЁ]+\b', text)

    # Step 3: Filter the words based on the conditions
    result_words = []
    for word in words:
        # Use the lowercase version for dictionary lookup
        word_lower = word.lower()

        if word_lower in word_data:
            syllables, stress_position = word_data[word_lower]

            # Condition 1: Exclude one-syllable words
            if syllables > 1:
                # Condition 2: Exclude words with the accent on the first syllable
                if stress_position != 1:
                    result_words.append(word)

    # Step 4: Print the final comma-separated list
    print(",".join(result_words))

solve_russian_stress_task()