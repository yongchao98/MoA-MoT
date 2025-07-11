import re

def solve_russian_stress_task():
    """
    Finds words in a given Russian text that are multi-syllable
    and do not have the stress on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # This dictionary acts as our linguistic knowledge base.
    # It maps lowercase words to the 0-indexed number of the stressed syllable.
    # For example, in 'шашлык' (шаш-лы́к), the stress is on the 2nd syllable, so its index is 1.
    stress_knowledge = {
        'шашлык': 1,    # шаш-лы́к
        'запах': 0,     # за́-пах
        'горелым': 1,   # го-ре́-лым
        'вскоре': 0,    # вско́-ре
        'прибежал': 2,  # при-бе-жа́л
        'сосед': 1,      # со-се́д
        'видно': 0,      # ви́д-но
        'дошёл': 1,     # до-шёл
        'его': 1,       # е-го́
        'дома': 0       # до́-ма (genitive case of 'дом', meaning 'of the house')
    }

    russian_vowels = "аеёиоуыэюя"

    # Use regex to find all sequences of Cyrillic letters (words)
    words = re.findall(r'[а-яА-ЯёЁ]+', text)

    result_words = []
    for word in words:
        # Work with a consistent lowercase version for lookups
        clean_word = word.lower()

        # 1. Filter out one-syllable words
        syllable_count = sum(1 for char in clean_word if char in russian_vowels)
        if syllable_count <= 1:
            continue

        # 2. Check if the stress is not on the first syllable
        if clean_word in stress_knowledge:
            stressed_syllable_index = stress_knowledge[clean_word]
            # The first syllable is index 0. We need words where the stress is on a later syllable.
            if stressed_syllable_index > 0:
                result_words.append(clean_word)

    print(", ".join(result_words))

solve_russian_stress_task()
<<<шашлык, горелым, прибежал, сосед, дошёл, его>>>