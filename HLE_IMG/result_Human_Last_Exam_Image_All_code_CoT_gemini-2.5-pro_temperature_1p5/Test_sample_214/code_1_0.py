def solve_poetic_form():
    """
    This function analyzes the syllable count of a given line of poetry
    to help identify the poetic form.
    """
    print("Analyzing the phrase 'nacreous wavers' to determine the poetic form.")
    print("The primary hypothesis is that the poem follows a syllabic form like a Haiku (5-7-5 syllables per line).\n")

    # A dictionary to store the syllable count for each word.
    syllable_counts = {
        "nacreous": 3,
        "wavers": 2
    }

    line_words = ["nacreous", "wavers"]
    total_syllables = 0
    calculation_string = []

    # Calculate the total syllables and build the output string.
    for word in line_words:
        count = syllable_counts.get(word, 0)
        total_syllables += count
        calculation_string.append(str(count))

    # Print the detailed breakdown of the calculation.
    print("Syllable calculation:")
    final_equation = f"{' + '.join(calculation_string)} = {total_syllables}"
    print(final_equation)

    # Print the conclusion based on the syllable count.
    print(f"\nThe line has {total_syllables} syllables.")
    print("This 5-syllable structure is the characteristic final line of a Haiku.")
    print("Therefore, the poetic form adhered to is the Haiku.")


solve_poetic_form()